### Run single task for electrolyte conductivity

import warnings
import os
from sklearn.model_selection import train_test_split
import configparser
import psycopg2
import argparse
from semiempirical import *
from utils import *

parser = configparser.ConfigParser()
parser.read("db.ini")
conn = psycopg2.connect(**parser["postgres"])

warnings.filterwarnings("ignore")

arg_parser = argparse.ArgumentParser(
    description="Electrolyte model parameters definition"
)

arg_parser.add_argument("--solvent1", default='EC', type=str, help="Input Solvent 1 name", required=True)
arg_parser.add_argument("--solvent2", default='PC', type=str, help="Input Solvent 2 name", required=True)
arg_parser.add_argument("--salt", default='LiAsF6', type=str, help="Input Salt name", required=True)
arg_parser.add_argument("--function_name", default='Landesfeind2019', choices=['Landesfeind2019', 'Weito2020', 'Kim2020'],
                        type=str, help="Input function name")

arg_parser.add_argument("--seed", default=100, type=int, help="Input seed for initial guess in curve fitting")
arg_parser.add_argument("--test_size", default=0.2, type=float, help="Input test size ratio for train test split")
arg_parser.add_argument("--random_split_seed", default=10, type=int, help="Input random seed for train test split")

args = arg_parser.parse_args()

## System definition
s1 = args.solvent1
s2 = args.solvent2
salt = args.salt
function_name = args.function_name

## Hyperparameter definition
seed = args.seed
test_size = args.test_size
random_split_seed = args.random_split_seed

## function info
function_info = conductivity_function_chosen(function_name, seed=seed)
curve_fit_function = function_info["curve_fit_function"]
initial_guess = function_info["initial_guess"]

## Data selection
df = system_chosen(s1, s2, salt)
X, Y = parse_dataset(df)

## Data splitting
xtrain, xtest, ytrain, ytest = train_test_split(X, Y, test_size=test_size, random_state=random_split_seed)

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

