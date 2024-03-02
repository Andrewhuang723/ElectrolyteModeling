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
engine = create_engine(f"postgresql://{parser['postgres']['user']}:{parser['postgres']['password']}@{parser['postgres']['dbname']}:5432/{parser['postgres']['dbname']}")

cols = ['index', 'Notes', 'Compound Notebook Name', 'Date Copied to Master Database', 'Data Recorded By', 'Seshadri Group', 'Oyaizu Group', 'DOI', 
'Solvent 1 (S1)', 'S1 SMILES', 'S1 SMILES 1', 'S1 SMILES 2', 'S1 SMILES 3', 'S1 Big SMILES', 'S1 End Group 1', 'S1 End Group 2', 'S1 Mol %', 'S1 Weight %', 'S1 SMILES 1 Mol %', 'S1 SMILES 2 Mol %', 'S1 SMILES 3 Mol %', 'S1 Molar Ratio', 'S1 Weight Ratio', 'S1 Mw', 'S1 Mn', 'S1 Mn or Mw', 'S1 Mw/Mn', 'S1 DOP', 'S1 Molar Mass (g/mol)', 'S1 Type', 'S1 Block Mol Weight', 
'Solvent 2 (S2)', 'S2 SMILES', 'S2 SMILES 1', 'S2 SMILES 2', 'S2 Mol %', 'S2 Weight %', 'S2 SMILES 1 Mol %', 'S2 SMILES 2 Mol %', 'S2 Molar Ratio', 'S2 Weight Ratio', 'S2 Mw', 'S2 Mn', 'S2 Mn or Mw', 'S2 Mw/Mn', 'S2 DOP', 'S2 Molar Mass (g/mol)', 'S2 Type', 
'Solvent 3 (S3)', 'S3 SMILES', 'S3 SMILES 1', 'S3 Mol %', 'S3 Weight %', 'S3 Molar Ratio', 'S3 Weight Ratio', 'S3 Mw', 'S3 Mn', 'S3 Mn or Mw', 'S3 Mw/Mn', 'S3 DOP', 'S3 Molar Mass (g/mol)', 'S3 Type', 
'Solvent 4 (S4)', 'S4 SMILES', 'S4 Mol %', 'S4 Weight %', 'S4 Molar Ratio', 'S4 Weight Ratio', 'S4 Mw', 'S4 Mn', 'S4 Mn or Mw', 'S4 Mw/Mn', 'S4 DOP', 'S4 Molar Mass (g/mol)', 'S4 Type', 
'Salt 1', 'Salt1 SMILES', 'Salt1 Mol %', 'Salt1 Weight %', 'Salt1 Molality (mol salt/kg polymer)', 'Salt1 Molar Mass (g/mol)',
'Salt 2', 'Salt2 SMILES', 'Salt2 Mol %', 'Salt2 Weight %', 'Salt2 Molality (mol salt/kg polymer)', 'Salt2 Molar Mass (g/mol)', 
'Inorganic Material 1 (IM1)', 'IM1 Weight %', 'IM1 Particle Size (nm)', 
'Tg (oC)', 'Heating Rate (K/min)', 'Temperature (oC)', 
'Conductivity (S/cm)', 'log Conductivity (S/cm)']

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

### SELECT database based on T8 with PEGDM-family polymer electrolyte
SELECT_T8_BASED_ELECTROLYTE = """
SELECT "index",
        "Solvent 1 (S1)", "S1 Weight %",
       "Solvent 2 (S2)", "S2 Weight %",
       "Salt 1",  "Salt1 Molality (mol salt/kg polymer)",
       "Temperature (oC)","Conductivity (S/cm)"
FROM polymerelectrolytedata
WHERE ("Solvent 1 (S1)" IN ('PC', 'EC', 'EMC', 'DMC')
AND "Solvent 2 (S2)" IN ('PC',  'EC', 'EMC', 'DMC'))
OR "Solvent 1 (S1)" LIKE 'PEGDM%'
AND "Salt 1" IN (
    SELECT "Salt 1"
    FROM polymerelectrolytedata
    GROUP BY "Salt 1"
    HAVING COUNT("Salt 1") > 50
);
"""


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
    AND "Salt 1" = '{salt}');
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
    AND "Salt 1" = '{salt}') 
OR (
    "Solvent 1 (S1)" = '{s2}'
    AND "Solvent 2 (S2)" = '{s1}'
    AND "Salt 1" = '{salt}' 
);
"""


SELECT_COLUMNS = """
SELECT column_name
FROM information_schema.columns
WHERE table_schema = 'public' 
AND table_name = 'polymerelectrolytedata';"""


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


def get_components(component_name: str, style="solvent"):
    with psycopg2.connect(**parser["postgres"]) as conn:
        with conn.cursor() as cur:
            if style == "solvent":
                cur.execute(SELECT_SOLVENT(solvent_name=component_name))
            elif style == "salt":
                cur.execute(SELECT_SALT(salt_name=component_name))
            else:
                raise Exception("No component style found")
            solvent_infos = cur.fetchall()
    return solvent_infos[0]

def get_index():
    with psycopg2.connect(**parser["postgres"]) as conn:
        with conn.cursor() as cur:
            cur.execute("""SELECT MAX("index") FROM polymerelectrolytedata""")
            max_idx = cur.fetchall()[0][0]
    return max_idx + 1


if __name__ == "__main__":

    ## Read original data into DB
    # ReadOriginalData2Postgres(url="https://raw.githubusercontent.com/learningmatter-mit/Chem-prop-pred/main/data/PolymerElectrolyteData.csv",
    #                             table_name="polymerelectrolytedata")
    PC_INFO = get_components("PC", style="solvent")
    EC_INFO = get_components("EC", style="solvent")
    EMC_INFO = get_components("EMC", style="solvent")
    DMC_INFO = get_components("DMC", style="solvent")
    LiPF6_INFO = get_components("LiPF6", style="salt")