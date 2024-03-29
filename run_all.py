### This is the run script for batch electrolyte systems runnning
### 2024/3/2 Andrew.Huang

import warnings
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.optimize import differential_evolution, curve_fit
import configparser
import psycopg2
from query import *
from semiempirical import *
from utils import *

parser = configparser.ConfigParser()
parser.read("db.ini")
conn = psycopg2.connect(**parser["postgres"])

## Update ethylene carbonate to EC & propylene carbonate to PC
conn.cursor().execute(UPDATE_SOLVENT_NAME)


### Select the unique combination of (solvent 1, solvent 2, salt)
SELECT_DISTINCT_S1_S2_SALT = """
SELECT COUNT(t1.system), t1.system
FROM (
SELECT CONCAT("Solvent 1 (S1)", '_',
       "Solvent 2 (S2)", '_',
       "Salt 1") AS system
FROM polymerelectrolytedata
WHERE ("Solvent 1 (S1)" IN ('PC', 'EC', 'EMC', 'DMC')
AND "Solvent 2 (S2)" IN ('PC',  'EC', 'EMC', 'DMC'))
OR "Solvent 1 (S1)" LIKE 'PEGDM%') t1
GROUP BY t1.system
HAVING COUNT(t1.system) > 10;
"""

## Select distinct sets of (solvent 1, solvent 2, salt)
electrolyte_systems = pd.read_sql_query(SELECT_DISTINCT_S1_S2_SALT, conn).values

electrolytes_sets = set()
warnings.filterwarnings("ignore")

total_results_df = None
counter = 0

for electrolyte_system in electrolyte_systems:
    ## [70, "S1-S2-Salt"]
    s1, s2, salt = electrolyte_system[-1].split("_")
    s1_dielectric_constant = dielectric_constant(s1)
    s2_dielectric_constant = dielectric_constant(s2)
    
    ## To check if the set exists or not
    if (s1, s2, salt) not in electrolytes_sets:
        electrolytes_sets.add((s1, s2, salt))
        electrolytes_sets.add((s2, s1, salt))
    else:
        continue
    
    print(f"\n{s1} || {s2} || {salt} :")

    try:
        ## chose the system where (Solvent 1 == s1 AND Solvent 2 == s2 AND Salt1 == salt)
        #                      or (Solvent 1 == s2 AND Solvent 2 == s1 AND Salt1 == salt)
        df = system_chosen(s1, s2, salt)
        
        # EC melting point ~ 30C
        if s1 == "EC" or s2 == "EC":
            df = df[df["Temperature (oC)"] > 34]

    except Exception as e:
        print(f"{e}")
        continue
    
    function_name = "Weito2020"
    # function_name = "Kim2011"
    # function_name = "Landesfeind2019"
    # function_name = "Ai2020"
    # function_name = "Ecker2015"
    # function_name = "Capiglia1999"
    # function_name = "Nyman2008_arrhenius"
    # function_name = "Prada2013"
    # function_name = "Ramadass2004"
    # function_name = "Valoen2005"

    seed = 100
    function_info = conductivity_function_chosen(function_name, seed=seed)
    curve_fit_function = function_info["curve_fit_function"]
    initial_guess = function_info["initial_guess"]

    X, Y = parse_dataset(df)
    
    try:
        xtrain, xtest, ytrain, ytest = train_test_split(X, Y, test_size=0.2, random_state=10)
    except ValueError:
        print("Not enough samples for train test split")
        continue
    # xtrain, xtest are tuple in the form of (solvent1, solvent2, salt, temp)

    try:
        parameter = run_fitting(xtrain, ytrain, 
                                curve_fit_function=curve_fit_function, 
                                initial_guess=initial_guess)
    except Exception as e:
        print(e)
        continue
    
    ypred = evaluate(xtest, ytest, 
            curve_fit_function=curve_fit_function, 
            parameter=parameter)
    
    ### Save results to csv
    results_df = xtest.copy()
    results_df["Conductivity_test"] = ytest
    results_df["Conductivity_pred"] = ypred
    results_df["error"] = np.abs(ypred - ytest) / ytest
    if not os.path.exists(f"./data/semiempirical_TEMP/{function_name}/{s1}_{s2}_{salt}"):
        os.makedirs(f"./data/semiempirical_TEMP/{function_name}/{s1}_{s2}_{salt}", exist_ok=True)
    results_df.to_csv(f"./data/semiempirical_TEMP/{function_name}/{s1}_{s2}_{salt}/results.csv", index=False)
    pd.DataFrame(parameter).to_csv(f"./data/semiempirical_TEMP/{function_name}/{s1}_{s2}_{salt}/parameter.csv", index=False)

    ### Plots each results
    ax = start_plot(style="darkgrid", dpi=180)
    R2_plot(reconstruct_y_pred=results_df["Conductivity_pred"].values, 
                reconstruct_y_test=results_df["Conductivity_test"].values, 
                error=results_df["error"].mean(),
                accuracy=len(results_df[results_df["error"] < 0.1]) / len(results_df),
                prop_name_T="Conductivity " + r"($\frac{mS}{cm}$)")
    plt.savefig(f"./data/semiempirical_TEMP/{function_name}/{s1}_{s2}_{salt}/R2_plot.png")

    if counter == 0:
        total_results_df = results_df.copy()
    else:
        total_results_df = pd.concat([total_results_df, results_df], axis=0)
    counter += 1

### Plots total results
total_results_df.to_csv(f"./data/semiempirical_TEMP/{function_name}/results.csv", index=False)
R2_plot(reconstruct_y_pred=total_results_df["Conductivity_pred"].values, 
            reconstruct_y_test=total_results_df["Conductivity_test"].values, 
            error=total_results_df["error"].mean(),
            accuracy=len(total_results_df[total_results_df["error"] < 0.1]) / len(total_results_df),
            prop_name_T="Conductivity " + r"($\frac{mS}{cm}$)")
plt.savefig(f"./data/semiempirical_TEMP/{function_name}/R2_plot.png")