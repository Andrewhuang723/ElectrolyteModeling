import pandas as pd
import numpy as np
import configparser
import psycopg2
import os
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from scipy.optimize import differential_evolution, curve_fit
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict


parser = configparser.ConfigParser()
parser.read("db.ini")
conn = psycopg2.connect(**parser["postgres"])


def conductivity_function_chosen(name: str, seed=42) -> Dict:
    """
    Returns the initial guess value for curve fit
    """
    np.random.seed(seed)
    if name == "Landesfeind2019":
        return {
            "curve_fit_function": curve_fit_function_Landesfeind2019,
            "function": Landesfeind2019_ElectrolyteConductivity,
            "initial_guess": np.random.rand(6)
        }
    elif name == "Weito2020":
        return {
            "curve_fit_function": curve_fit_function_Weito2020,
            "function": Weito2020_ElectolyteConductivity,
            "initial_guess": np.random.rand(5)
        }

    elif name == "Kim2011":
        return {
            "curve_fit_function":  curve_fit_function_Kim2011,
            "function": Kim2011_ElectrolyteConductivity,
            "initial_guess": np.random.rand(3)
        }
    else:
        raise Exception(f"Function {name} is not found.\n Please input 'Landesfeind2019' or Weito2020' ")


def parse_dataset(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convert original dataframe
    Conductivity [S/cm] -> [ms/cm]
    Temperature [oC] -> [K]
    and return X & Y
    """
    X = df.drop(columns=["index", "Conductivity (S/cm)"], axis=1)
    X["Temperature (K)"] = X["Temperature (oC)"].apply(lambda t: t + 273.15)
    X.drop(columns=["Temperature (oC)"], inplace=True)
    Y = df["Conductivity (S/cm)"] * 1000
    return X, Y


def dielectric_constant(name):
    if name == "EC":
        return 90.36
    elif name == "EMC":
        return 2.40
    elif name == "DMC":
        return 3.108
    elif name == "PC":
        return 64.95
    elif name == "PEGDM":
        return 12.5

def Weito2020_ElectolyteConductivity(solvent_1, solvent_2, salt, temp, p1, p2, p3, p4, n):
    return (p1 * temp + p2) * np.power(salt, n) * np.exp((-p3 * salt)/(temp - p4))

def Kim2011_ElectrolyteConductivity(solvent_1, solvent_2, salt, temp, p1, p2, p3):
    return p1 * np.exp(-798 / temp) * (salt) ** 3 + p2 * np.exp(-1080 / temp) * (salt) ** 2 + p3 * np.exp(-1440 / temp) * (salt)

def Landesfeind2019_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5, p6):
    c = c_e  # mol.m-3 -> mol.l
    # p1, p2, p3, p4, p5, p6 = coeffs
    A = p1 * (1 + (T - p2))
    B = 1 + p3 * np.sqrt(c) + p4 * (1 + p5 * np.exp(1000 / T)) * c
    C = 1 + c ** 4 * (p6 * np.exp(1000 / T))
    sigma_e = A * c * B / C  # mS.cm-1

    return sigma_e / 10


def curve_fit_function_Weito2020(X: tuple, *args):
    return Weito2020_ElectolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Landesfeind2019(X: tuple, *args):
    return Landesfeind2019_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Kim2011(X: tuple, *args):
    return Kim2011_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)


def parse_XY(X: pd.DataFrame, Y: pd.DataFrame, s2=None) -> Tuple[Tuple[np.array, np.array, np.array, np.array], np.array]:
    """
    Return the values 
    Xs: tuple[], prepared for scipy.optimize.curve_fit
    Ys: 
    """
    solvent_1_value = X["S1 Weight %"].values
    if not s2:
        solvent_2_value = np.zeros(len(solvent_1_value))
    else:
        solvent_2_value = X["S2 Weight %"].values
    salt_value = X["Salt1 Molality (mol salt/kg polymer)"].values
    temperature = X["Temperature (K)"].values


    Xs = (solvent_1_value, solvent_2_value, salt_value, temperature)

    Ys = Y.values
    return Xs, Ys


def run_fitting(xtrain: pd.DataFrame, ytrain: pd.DataFrame, curve_fit_function: callable, initial_guess):
    xtrain_prepared, ytrain_prepared = parse_XY(xtrain, ytrain)

    n_features = len(xtrain_prepared)
    n_samples = xtrain_prepared[0].shape[0]
    if n_features > n_samples:
        raise Exception(f"n_features: {n_features} > n_samples: {n_samples}")

    p0=initial_guess
    try:
        parameter, covariace = curve_fit(curve_fit_function, xtrain_prepared, ytrain_prepared, p0=initial_guess)
        return parameter
    except RuntimeError:
        raise Exception("curve_fit encounter RuntimeError")


def evaluate(xtest: pd.DataFrame, ytest:pd.DataFrame, curve_fit_function: callable, parameter: np.array):
    xtest_prepared, ytest_prepared = parse_XY(xtest, ytest)
    ypred = curve_fit_function(xtest_prepared, *parameter)
    r2 = r2_score(ypred, ytest_prepared)
    error = np.abs(ypred - ytest_prepared) / ytest_prepared
    accuracy = len(np.where(error < 0.1)[0]) / len(ytest_prepared)

    print(r2, accuracy)
    return ypred


def start_plot(figsize=(10, 8), style = 'whitegrid', dpi=100):
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(1,1)
    plt.tight_layout()
    with sns.axes_style(style):
        ax = fig.add_subplot(gs[0,0])
    return ax

def R2_plot(reconstruct_y_pred, reconstruct_y_test, error, accuracy, prop_name_T, ax=None):
    if ax is None:
        ax = start_plot(style='darkgrid', dpi=180)
    
    sns.set(rc={'font.family': 'Times New Roman'})
    ax.scatter(reconstruct_y_test.reshape(-1), reconstruct_y_pred.reshape(-1), color='darkorange', edgecolor='navy',
               label=r'$R^2$' + ": %.4f" % r2_score(y_true=reconstruct_y_test, y_pred=reconstruct_y_pred) + '\n' +
                     "MAPE" + ": %.4f" % error + '\n' +
                     "Accuracy" + ": %.4f" % accuracy)
    ymin = min(np.min(reconstruct_y_test), np.min(reconstruct_y_pred)) - 0.1
    ymax = max(np.max(reconstruct_y_test), np.max(reconstruct_y_pred)) + 0.1
    lim = [ymin, ymax]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, c='brown', ls='--', label=r'$y=\hat y, $' + 'identity')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=15)
    plt.xlabel('TRUE %s' % prop_name_T, fontsize=20, font="Times New Roman")
    plt.ylabel('PREDICTED %s' % prop_name_T, fontsize=20, font="Times New Roman")
    plt.xticks(fontsize=20, font="Times New Roman")
    plt.yticks(fontsize=20, font="Times New Roman")
    plt.title('%s Testing: %d' % (prop_name_T, len(reconstruct_y_test)), fontsize=20, font="Times New Roman")
    return ax


if __name__ == "__main__":

    ### Run combination
    import warnings
    from query import UPDATE_SOLVENT_NAME, SELECT_DISTINCT_S1_S2_SALT, SELECT_BINARY_SYSTEMS, SELECT_SINGLE_SYSTEMS, system_chosen
    
    ## Update ethylene carbonate to EC & propylene carbonate to PC
    conn.cursor().execute(UPDATE_SOLVENT_NAME)
    
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
        except Exception as e:
            print(f"{e}")
            continue
        
        # function_name = "Weito2020"
        # function_name = "Kim2011"
        function_name = "Landesfeind2019"
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
    
        if counter == 0:
            total_results_df = results_df.copy()
        else:
            total_results_df = pd.concat([total_results_df, results_df], axis=0)
        counter += 1

    ### Plots total results
    total_results_df.to_csv(f"./data/semiempirical/{function_name}/results.csv", index=False)
    R2_plot(reconstruct_y_pred=total_results_df["Conductivity_pred"].values, 
                reconstruct_y_test=total_results_df["Conductivity_test"].values, 
                error=total_results_df["error"].mean(),
                accuracy=len(total_results_df[total_results_df["error"] < 0.1]) / len(total_results_df),
                prop_name_T="Conductivity " + r"($\frac{mS}{cm}$)")
    plt.savefig(f"./data/semiempirical/{function_name}/R2_plot.png")
