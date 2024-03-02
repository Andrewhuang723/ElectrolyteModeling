import warnings
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.optimize import differential_evolution, curve_fit
import configparser
import psycopg2
from semiempirical import *
from utils import *

parser = configparser.ConfigParser()
parser.read("db.ini")
conn = psycopg2.connect(**parser["postgres"])


warnings.filterwarnings("ignore")


## System definition
s1 = "EC" # solvent 1
s2 = "PC" # solvent 2
salt = "LiAsF6" # salt

## Function defintion
function_name = "Landesfeind2019"

## Hyperparameter definition
seed = 100
function_info = conductivity_function_chosen(function_name, seed=seed)
curve_fit_function = function_info["curve_fit_function"]
initial_guess = function_info["initial_guess"]

## Data selection
df = system_chosen(s1, s2, salt)
X, Y = parse_dataset(df)

## Data splitting
xtrain, xtest, ytrain, ytest = train_test_split(X, Y, test_size=0.2, random_state=10)

## Run fitting -> return parameter for the fitted function
parameter = run_fitting(xtrain, ytrain, 
                        curve_fit_function=curve_fit_function, 
                        initial_guess=initial_guess)

## Evaluation
ypred = evaluate(xtest, ytest, 
        curve_fit_function=curve_fit_function, 
        parameter=parameter)


### Save results to csv
results_df = xtest.copy()
results_df["Conductivity_test"] = ytest
results_df["Conductivity_pred"] = ypred
results_df["error"] = np.abs(ypred - ytest) / ytest
if not os.path.exists(f"./data/semiempirical/{function_name}/{s1}_{s2}_{salt}"):
    os.mkdir(f"./data/semiempirical/{function_name}/{s1}_{s2}_{salt}")
results_df.to_csv(f"./data/semiempirical/{function_name}/{s1}_{s2}_{salt}/results.csv", index=False)
pd.DataFrame(parameter).to_csv(f"./data/semiempirical/{function_name}/{s1}_{s2}_{salt}/parameter.csv", index=False)

### Plots each results
ax = start_plot(style="darkgrid", dpi=180)
R2_plot(reconstruct_y_pred=results_df["Conductivity_pred"].values, 
            reconstruct_y_test=results_df["Conductivity_test"].values, 
            error=results_df["error"].mean(),
            accuracy=len(results_df[results_df["error"] < 0.1]) / len(results_df),
            prop_name_T="Conductivity " + r"($\frac{mS}{cm}$)")
plt.savefig(f"./data/semiempirical/{function_name}/{s1}_{s2}_{salt}/R2_plot.png")