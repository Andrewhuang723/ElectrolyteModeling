import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import configparser
import psycopg2
from typing import Tuple, Dict

parser = configparser.ConfigParser()
parser.read("db.ini")
conn = psycopg2.connect(**parser["postgres"])


#constants
R = 8.314


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
    elif name == "Ai2020":
        return {
            "curve_fit_function":  curve_fit_function_Ai2020,
            "function": Ai2020_ElectrolyteConductivity,
            "initial_guess": np.random.rand(8)
        }
    elif name == "Ecker2015":
        return {
            "curve_fit_function":  curve_fit_function_Ecker2015,
            "function": Ecker2015_ElectrolyteConductivity,
            "initial_guess": np.random.rand(5)
        }
    elif name == "Capiglia1999":
        return {
            "curve_fit_function":  curve_fit_function_Capiglia1999,
            "function": Capiglia1999_ElectrolyteConductivity,
            "initial_guess": np.random.rand(5)
        }
    elif name == "Nyman2008_arrhenius":
        return {
            "curve_fit_function":  curve_fit_function_Nyman2008_arrhenius,
            "function": Nyman2008_arrhenius_ElectrolyteConductivity,
            "initial_guess": np.random.rand(4)
        }
    elif name == "Prada2013":
        return {
            "curve_fit_function":  curve_fit_function_Prada2013,
            "function": Prada2013_ElectrolyteConductivity,
            "initial_guess": np.random.rand(6)
        }
    elif name == "Ramadass2004":
        return {
            "curve_fit_function":  curve_fit_function_Ramadass2004,
            "function": Ramadass2004_ElectrolyteConductivity,
            "initial_guess": np.random.rand(6)
        }
    elif name == "Valoen2005":
        return {
            "curve_fit_function":  curve_fit_function_Valoen2005,
            "function": Valoen2005_ElectrolyteConductivity,
            "initial_guess": np.random.rand(8)
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

def Ai2020_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5, p6, p7, p8):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration.
    Concentration should be in dm3 in the function.

    Sourced from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Ai2020.py
    """

    sigma_e = (
        c_e
        * (
            (p1 + p2 * c_e + p3 * c_e**2)
            + (p4 + p5 * c_e + p6 * c_e**2) * T
            + (p8 * + p7 * c_e) * T**2
        )
        ** 2
    )

    return sigma_e


def Ecker2015_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration [1, 2, 3].

    Sourced from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Ecker2015.py
    """


    # value at T = 296K
    sigma_e_296 = p1 * c_e **3 - p2 * c_e **2 + p3 * c_e + p4

    # add temperature dependence
    # E_k_e = 1.71e4
    E_k_e = p5
    C = 296 * np.exp(E_k_e / (R * 296))
    sigma_e = C * sigma_e_296 * np.exp(-E_k_e / (R * T)) / T

    return sigma_e


def Capiglia1999_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration. The original
    data is from [1]. The fit is from Dualfoil [2].

    Sourced from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Marquis2019.py
    """

    sigma_e = (
        p1
        + p2 * (c_e)
        - p3 * (c_e) ** 2
        + p4 * (c_e) ** 3
    )

    E_k_e = p5
    arrhenius = np.exp(E_k_e / R * (1 / 298.15 - 1 / T))

    return sigma_e * arrhenius



def Nyman2008_arrhenius_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4):
    """
    Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1], with Arrhenius temperature dependence added from [2].

    Source from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/OKane2022.py
    """

    sigma_e = (
        p1 * (c_e) ** 3 - p2 * (c_e) ** 1.5 + p3 * (c_e)
    )

    # Nyman et al. (2008) does not provide temperature dependence
    # So use temperature dependence from Ecker et al. (2015) instead

    E_sigma_e = p4
    arrhenius = np.exp(E_sigma_e / R * (1 / 298.15 - 1 / T))

    return sigma_e * arrhenius


def Prada2013_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5, p6):
    """
    Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from :footcite:`Prada2013`.

    Source from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Prada2013.py
    """
    # convert c_e from mol/m3 to mol/L
    # c_e = c_e / 1e6

    sigma_e = (
        p1
        + p2 * c_e
        + p3 * c_e**2
        + p4 * c_e**3
        + p5 * c_e**4
    )

    # E_sigma_e = 17000
    # arrhenius = np.exp(E_sigma_e / R * (1 / 298.15 - 1 / T))
    E_k_e = p6
    arrhenius = np.exp(E_k_e / R * (1 / 298.15 - 1 / T))

    return sigma_e * arrhenius

def Ramadass2004_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5, p6):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration.
    Concentration should be in dm3 in the function.

    Source from PyBamm

    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Ramadass2004.py
    """
    # mol.m-3 to mol.dm-3, original function is likely in mS/cm
    # The function is not in Arora 2000 as reported in Ramadass 2004

    # cm = 1e-6 * c_e  # here it should be only 1e-3

    sigma_e = (
        p1
        + p2 * c_e
        + p3 * (c_e**2)
        + p4 * (c_e**3)
        + p5 * (c_e**8) * (c_e**4)
    )  # and here there should not be an exponent

    E_k_e = p6
    arrhenius = np.exp(E_k_e / R * (1 / 298.15 - 1 / T))

    return sigma_e * arrhenius


def Valoen2005_ElectrolyteConductivity(solvent_1, solvent_2, c_e, T, p1, p2, p3, p4, p5, p6, p7, p8):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration, from [1]
    (eqn 17)

    Source from PyBamm
    https://github.com/pybamm-team/PyBaMM/blob/bc31cacb5e7fcd2d1d16dd6d58d63e2db4474a7b/pybamm/input/parameters/lithium_ion/Xu2019.py
    """
    # mol/m3 to molar
    # c_e = c_e / 1000
    # mS/cm to S/m
    return (
        c_e
        * (
            (p1 + p2 * T + p3 * T**2)
            + c_e * (p4 + p5 * T + p6 * T**2)
            + c_e**2 * (p7 + p8 * T)
        )
        ** 2
    )


def curve_fit_function_Weito2020(X: tuple, *args):
    return Weito2020_ElectolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Landesfeind2019(X: tuple, *args):
    return Landesfeind2019_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Kim2011(X: tuple, *args):
    return Kim2011_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Ai2020(X: tuple, *args):
    return Ai2020_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Ecker2015(X: tuple, *args):
    return Ecker2015_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Capiglia1999(X: tuple, *args):
    return Capiglia1999_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Nyman2008_arrhenius(X: tuple, *args):
    return Nyman2008_arrhenius_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Prada2013(X: tuple, *args):
    return Prada2013_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Ramadass2004(X: tuple, *args):
    return Ramadass2004_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)

def curve_fit_function_Valoen2005(X: tuple, *args):
    return Valoen2005_ElectrolyteConductivity(X[0], X[1], X[2], X[3], *args)




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

    print("R2: %.4f" % r2, "Accuracy: %.4f" % accuracy)
    return ypred