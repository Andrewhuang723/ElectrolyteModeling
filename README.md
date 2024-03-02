# ElectrolyteModeling
Build model for electrolyte system in Lithium Ion Battery

---
## Semiempirical Model for Li+ conductivity of electrolytes

### Database

- Data source: [Chem-prop-pred](https://github.com/learningmatter-mit/Chem-prop-pred/blob/main/data/PolymerElectrolyteData.csv)[1]

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
