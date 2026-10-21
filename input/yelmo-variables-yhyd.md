# yhyd

| id | variable          | dimensions        | units        | long_name                                                | prec  |
|----|-------------------|-------------------|--------------|----------------------------------------------------------|-------|
|  1 | hyd_W_til         | xc, yc            | m            | Till water storage thickness (bucket)                    | 1e-3  |
|  2 | hyd_W_til_max     | xc, yc            | m            | Per-cell cap on till water storage                       | 1e-3  |
|  3 | hyd_dW_til_dt     | xc, yc            | m/s          | Rate of change of till water storage                     | 1e-4r |
|  4 | hyd_overflow      | xc, yc            | m/s          | Till-saturation spill rate (K24 source)                  | 1e-4r |
|  5 | hyd_W             | xc, yc            | m            | K24 distributed sheet thickness                          | 1e-4r |
|  6 | hyd_p_w           | xc, yc            | Pa           | Basal water pressure                                     | 1e-4r |
|  7 | hyd_q_x           | xc, yc            | m^2/s        | Basal water flux, x-component                            | 1e-4r |
|  8 | hyd_q_y           | xc, yc            | m^2/s        | Basal water flux, y-component                            | 1e-4r |
|  9 | hyd_N             | xc, yc            | Pa           | Effective pressure at the bed                            | 1e-4r |
| 10 | hyd_kappa         | xc, yc            | -            | K24 hydraulic transmissivity field                       | 1e-4r |
| 11 | hyd_q             | xc, yc            | m^2/s        | K24 distributed flux magnitude (Picard warm-start state) | 1e-4r |
| 12 | hyd_C_frz         | xc, yc            | m/s          | Freeze-on capacity at the bed (ice equivalent)           | 1e-4r |
| 13 | hyd_Q_diss        | xc, yc            | W/m^2        | Heat dissipated by basal water flow                      | 1e-4r |
| 14 | hyd_Q_sens        | xc, yc            | W/m^2        | Sensible heat of basal water flow                        | 1e-4r |
