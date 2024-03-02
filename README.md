# ElectrolyteModeling
Build model for electrolyte system in Lithium Ion Battery

---
## Semiempirical Model for Li+ conductivity of electrolytes

### Database

- Data source: [PolymerElectrolyteData](https://github.com/learningmatter-mit/Chem-prop-pred/blob/main/data/PolymerElectrolyteData.csv)[1]

- Queries: 
    
    We perform postgresql queries to request data from original database in `query.py`

    ```
    SELECT "Solvent 1 (S1)", "Salt 1", "Temperature (oC)", "Conductivity (S/cm)" FROM polymerelectrolytedata 
    WHERE "Solvent 1 (S1)" = 'EC'
    LIMIT 3;
    ```

    Return the value as following:

    |  Solvent 1 (S1) | Salt 1          | Temperature (oC) | Conductivity (S/cm) |
    | --------------- | --------------- | ---------------- | ------------------- |
    | EC              | LiTFSI          | 49.93288271      | 0.000281838         |
    | EC              | LITFSI          | 89.87674273      | 0.000432229         |
    | EC              | LITFSI          | 0.092633174	   | 0.0000941           |



### Models

- Semiempirical Model for elecrolyte conductivity ($\sigma_e$):

    $\sigma_e = f(m, T)$

    $m$: salt concentration (mol / kg polymer or M)

    $T$: temperature (K)

- We implemented 3 semiempirical models for electrolyte conductivity in `semiempirical.py`:

    1. `Weitoa2020` [2]
    2. `LandesFiend2019` [3]
    3. `Kim2011` [4]


### Run scripts

#### Arguments

- `solvent1`: `str` (Required)

    - Define the first solvent for the electrolyte system.

    - Default: `"EC"`

- `solvent2`: `str` (Required)
    
    - Define the second solvent for the electrolyte system.

    - Default: `"PC"`

- `salt`: `str` (Required)
    
    - Define the lithium salt for the electrolyte system.

    - Default: `"LiAsF6"`


- `function_name`: `str` (Optional)

    - Define the semiempircal model for the electrolyte conductivtiy model.

    - Default: `"Landefiend2019"`

- `seed`: `int` (Optional)

    - Define the seed of the initial guess for curve fitting

    - Default: `100`

- `test_size`: `float` (Optional)

    - Define the test size for `train_test_split` function in `scikit-learn`

    - Default: `0.2`

- `random_split_seed`: `int` (Optional)

    - Define the seed for `train_test_split` function in `scikit-learn`.

    - Default: `10`


Example command:

- `"Landesfeind2019"`
```
$ python3 run.py --solvent1 "EC" --solvent2 "PC" --salt "LiAsF6" --function_name "Landesfeind2019"
```

- `"Weito2020"`
```
$ python3 run.py --solvent1 "EC" --solvent2 "PC" --salt "LiAsF6" --function_name "Weito2020"
```

- `"Kim2011"`
```
$ python3 run.py --solvent1 "DMC" --solvent2 "EC" --salt "LiTDI" --function_name "Kim2011"
```

Outputs:

```
R2: 0.9700 Accuracy: 0.3684
```

Also, return the result files in `./data/semiempirical/Landesfeind2019/EC_PC_LiAsF6`

```
├── R2_plot.png
├── parameter.csv
└── results.csv
```


### Dashboard

- Visualization results and metrics in `semiempirical_dashboard.py`


---
[1] https://github.com/learningmatter-mit/Chem-prop-pred/tree/main

[2] Weitao Zhang, Xia Chen, Yan Wang, Lianying Wu, and Yangdong Hu
ACS Omega 2020 5 (35), 22465-22474
DOI: 10.1021/acsomega.0c03013

[3] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

[4] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.
