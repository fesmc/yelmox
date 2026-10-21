# ydata

| id | variable          | dimensions          | units       | long_name                                     | prec       |
|----|-------------------|---------------------|-------------|-----------------------------------------------|------------|
|  1 | pd_H_ice          | xc, yc              | m           | PD ice thickness                              | 1e-3       |
|  2 | pd_z_srf          | xc, yc              | m           | PD surface elevation                          | 1e-3       |
|  3 | pd_z_bed          | xc, yc              | m           | PD bedrock elevation                          | 1e-3       |
|  4 | pd_H_grnd         | xc, yc              | m           | PD overburden ice thickness                   | 1e-3       |
|  5 | pd_mask_bed       | xc, yc              | -           | PD mask                                       | -          |
|  6 | pd_ux_s           | xc, yc              | m/yr        | PD surface velocity in the x-direction        | 1e-3,1e-4r |
|  7 | pd_uy_s           | xc, yc              | m/yr        | PD surface velocity in the y-direction        | 1e-3,1e-4r |
|  8 | pd_uxy_s          | xc, yc              | m/yr        | PD surface velocity magnitude                 | 1e-3,1e-4r |
|  9 | pd_T_srf          | xc, yc              | K           | PD surface temperature                        | 1e-3       |
| 10 | pd_smb_ref        | xc, yc              | m/yr        | PD surface mass balance                       | 1e-4       |
| 11 | pd_depth_iso      | xc, yc, pd_time_iso | m           | PD depth of specific isochrones               | 1e-3       |
| 12 | pd_err_H_ice      | xc, yc              | m           | PD error in ice thickness                     | 1e-3       |
| 13 | pd_err_z_srf      | xc, yc              | m           | PD error in surface elevation                 | 1e-3       |
| 14 | pd_err_z_bed      | xc, yc              | m           | PD error in bedrock elevation                 | 1e-3       |
| 15 | pd_err_smb_ref    | xc, yc              | m/yr        | PD error in surface mass balance              | 1e-4       |
| 16 | pd_err_uxy_s      | xc, yc              | m/yr        | PD error in surface velocity magnitude        | 1e-3,1e-4r |
| 17 | pd_err_depth_iso  | xc, yc, pd_time_iso | m           | PD error in isochrone depth                   | 1e-3       |
