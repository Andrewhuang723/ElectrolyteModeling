import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
import configparser
import psycopg2
import os
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
        raise Exception(f"Function {name} is not found.\n Please input 'Landesfeind2019' or Weito2020' or 'Kim2011'")


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