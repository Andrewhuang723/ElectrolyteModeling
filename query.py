### Numerous of SQL query to request data from original database

import configparser
import pandas as pd
import psycopg2
import os
from sqlalchemy import create_engine
import numpy as np

parser = configparser.ConfigParser()

## Please create your own db.ini file
parser.read("db.ini")

## Create connection and engine
conn = psycopg2.connect(**parser["postgres"])
engine = create_engine(f"postgresql://{parser['postgres']['user']}:{parser['postgres']['password']}@{parser['postgres']['host']}:5432/{parser['postgres']['dbname']}")

### Update trivial solvent name to abbreviation
UPDATE_SOLVENT_NAME = """
UPDATE polymerelectrolytedata
SET "Solvent 1 (S1)" = CASE 
                            WHEN "Solvent 1 (S1)" = 'propylene carbonate' THEN 'PC'
                            WHEN "Solvent 1 (S1)" = 'ethylene carbonate' THEN 'EC'
                            ELSE "Solvent 1 (S1)"
                        END,
    "Solvent 2 (S2)" = CASE 
                            WHEN "Solvent 2 (S2)" = 'propylene carbonate' THEN 'PC'
                            WHEN "Solvent 2 (S2)" = 'ethylene carbonate' THEN 'EC'
                            ELSE "Solvent 2 (S2)"
                        END;
"""


### Select electrolyte system with only 1 solvent 
SELECT_SINGLE_SYSTEMS = lambda s1, salt: f"""
SELECT "index", 
    "Solvent 1 (S1)", "S1 Weight %",
       "Solvent 2 (S2)", "S2 Weight %",
       "Salt 1", "Salt1 Molality (mol salt/kg polymer)",
       "Temperature (oC)","Conductivity (S/cm)"
FROM polymerelectrolytedata
WHERE (
    "Solvent 1 (S1)" = '{s1}'
    AND "Salt 1" = '{salt}'
);
"""

## Select electrolyte system with binary solvents 
SELECT_BINARY_SYSTEMS = lambda s1, s2, salt: f"""
SELECT "index",
        "Solvent 1 (S1)", "S1 Weight %",
       "Solvent 2 (S2)", "S2 Weight %",
       "Salt 1", "Salt1 Molality (mol salt/kg polymer)",
       "Temperature (oC)","Conductivity (S/cm)"
FROM polymerelectrolytedata
WHERE (
    "Solvent 1 (S1)" = '{s1}'
    AND "Solvent 2 (S2)" = '{s2}'
    AND "Salt 1" = '{salt}'
) 
OR (
    "Solvent 1 (S1)" = '{s2}'
    AND "Solvent 2 (S2)" = '{s1}'
    AND "Salt 1" = '{salt}' 
);
"""


SELECT_SOLVENT = lambda solvent_name: f"""
SELECT "Solvent 1 (S1)", "S1 SMILES", "S1 Mn or Mw", 
"S1 Molar Mass (g/mol)"
FROM polymerelectrolytedata
WHERE "Solvent 1 (S1)" = '{solvent_name}' LIMIT 1;
"""

SELECT_SALT = lambda salt_name: f"""
SELECT "Salt 1", "Salt1 SMILES", "Salt1 Molar Mass (g/mol)"
FROM polymerelectrolytedata
WHERE "Salt 1" = '{salt_name}'
"""


def ReadOriginalData2Postgres(url, table_name):
    df = pd.read_csv(url)
    df.to_sql(table_name, engine, index=True, if_exists='replace')
    print(f"{url} is converted to SQL table {table_name}")


def system_chosen(s1, s2, salt):
    """
    System of electrolyte is composed of 3 components: solvent 1 & solvent 2 & salt
    If solvent 1 & solvent 2 were switched, then swap
    """
    if s2 is None or s2 == '':
        df = pd.read_sql_query(SELECT_SINGLE_SYSTEMS(s1, salt), con=conn)
    else:
        df = pd.read_sql_query(SELECT_BINARY_SYSTEMS(s1, s2, salt), con=conn)

        ## Swap values if S1 == S2 & S2 == S1

        condition = df["Solvent 1 (S1)"] == s2
        if condition.sum(axis=0) == 0:
            return df
        else:
            df.loc[condition, "S2 Weight %"], df.loc[condition, "S1 Weight %"] = df.loc[condition, "S1 Weight %"], df.loc[condition, "S2 Weight %"]
            df.loc[condition, "Solvent 2 (S2)"], df.loc[condition, "Solvent 1 (S1)"] = df.loc[condition, "Solvent 1 (S1)"], df.loc[condition, "Solvent 2 (S2)"]
    return df



if __name__ == "__main__":

    ## Read original data into DB
    # ReadOriginalData2Postgres(url="https://raw.githubusercontent.com/learningmatter-mit/Chem-prop-pred/main/data/PolymerElectrolyteData.csv",
    #                             table_name="polymerelectrolytedata")

    select = """
    SELECT *
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
    df = pd.read_sql_query(select, con=conn)
    print(df)