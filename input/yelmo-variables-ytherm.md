# ytherm

| id | variable          | dimensions        | units        | long_name                                     |
|----|-------------------|-------------------|--------------|-----------------------------------------------|
|  1 | enth              | xc, yc, zeta      | J kg^-1      | Ice enthalpy                                  |
|  2 | T_ice             | xc, yc, zeta      | K            | Ice temperature                               |
|  3 | omega             | xc, yc, zeta      | -            | Ice water content                             |
|  4 | T_pmp             | xc, yc, zeta      | K            | Pressure-corrected melting point              |
|  5 | T_prime           | xc, yc, zeta      | deg C        | Homologous ice temperature                    |
|  6 | f_pmp             | xc, yc            | -            | Fraction of cell at pressure melting point    |
|  7 | bmb_grnd          | xc, yc            | m/yr         | Grounded basal mass balance                   |
|  8 | Q_strn            | xc, yc, zeta      | J yr^-1 m^-3 | Internal strain heat production               |
|  9 | dQsdT             | xc, yc, zeta      | yr^-1        | Change of internal heat production w.r.t. T   |
| 10 | Q_b               | xc, yc            | mW m^-2      | Basal friction heat production                |
| 11 | Q_ice_b           | xc, yc            | mW m^-2      | Basal ice heat flux                           |
| 12 | T_prime_b         | xc, yc            | K            | Homologous temperature at the base            |
| 13 | cp                | xc, yc, zeta      | J kg^-1 K^-1 | Specific heat capacity                        |
| 14 | kt                | xc, yc, zeta      | J yr^-1 m^-1 K^-1 | Heat conductivity                             |
| 15 | H_cts             | xc, yc            | m            | Height of the CTS (cold-temperate surface)    |
| 16 | advecxy           | xc, yc, zeta      | J kg^-1 yr^-1 | Horizontal advection of enth (of T_ice in K yr^-1 if method=temp) |
| 17 | Q_rock            | xc, yc            | mW m^-2      | Heat flux from bedrock                        |
| 18 | T_rock            | xc, yc, zeta_rock | K            | Bedrock temperature                           |
| 19 | bmb_grnd_star     | xc, yc            | m/yr         | Grounded bmb of a base held at T_pmp (capacity rule) |
| 20 | bc_b              | xc, yc            | -            | Basal BC used: 0 none, 1 held at T_pmp, 2 flux |
| 21 | bmb_clamp         | xc, yc            | m/yr         | Freeze-on removed by the capacity clamp       |
| 22 | melt_int          | xc, yc            | m/yr         | Englacial water drained to the bed            |
