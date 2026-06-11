# yhyd

| id | variable          | dimensions        | units        | long_name                                     |
|----|-------------------|-------------------|--------------|-----------------------------------------------|
|  1 | hyd_W_til         | xc, yc            | m            | Till water storage thickness (bucket)         |
|  2 | hyd_W_til_max     | xc, yc            | m            | Per-cell cap on till water storage            |
|  3 | hyd_dW_til_dt     | xc, yc            | m/s          | Rate of change of till water storage          |
|  4 | hyd_overflow      | xc, yc            | m/s          | Till-saturation spill rate (K24 source)       |
|  5 | hyd_W             | xc, yc            | m            | K24 distributed sheet thickness               |
|  6 | hyd_p_w           | xc, yc            | Pa           | Basal water pressure                          |
|  7 | hyd_q_x           | xc, yc            | m^2/s        | Basal water flux, x-component                 |
|  8 | hyd_q_y           | xc, yc            | m^2/s        | Basal water flux, y-component                 |
|  9 | hyd_N             | xc, yc            | Pa           | Effective pressure at the bed                 |
| 10 | hyd_kappa         | xc, yc            | -            | K24 hydraulic transmissivity field            |
