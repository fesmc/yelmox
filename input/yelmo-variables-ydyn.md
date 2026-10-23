# ydyn

| id | variable          | dimensions       | units       | long_name                                     | prec       |
|----|-------------------|------------------|-------------|-----------------------------------------------|------------|
|  1 | ux                | xc, yc, zeta     | m/yr        | x-velocity                                    | 1e-3,1e-4r |
|  2 | uy                | xc, yc, zeta     | m/yr        | y-velocity                                    | 1e-3,1e-4r |
|  3 | uxy               | xc, yc, zeta     | m/yr        | Horizonal velocity magnitude                  | 1e-3,1e-4r |
|  4 | uz                | xc, yc, zeta_ac  | m/yr        | z-component velocity                          | 1e-3,1e-4r |
|  5 | uz_star           | xc, yc, zeta_ac  | m/yr        | z-velocity with corr. for thermal advection   | 1e-3,1e-4r |
|  6 | ux_bar            | xc, yc           | m/yr        | Depth-averaged x-velocity                     | 1e-3,1e-4r |
|  7 | uy_bar            | xc, yc           | m/yr        | Depth-averaged y-velocity                     | 1e-3,1e-4r |
|  8 | uxy_bar           | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  | 1e-3,1e-4r |
|  9 | ux_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged x-velocity            | 1e-3,1e-4r |
| 10 | uy_bar_prev       | xc, yc           | m/yr        | Previous depth-averaged y-velocity            | 1e-3,1e-4r |
| 11 | ux_b              | xc, yc           | m/yr        | Basal x-velocity                              | 1e-3,1e-4r |
| 12 | uy_b              | xc, yc           | m/yr        | Basal y-velocity                              | 1e-3,1e-4r |
| 13 | uz_b              | xc, yc           | m/yr        | Basal z-velocity                              | 1e-3,1e-4r |
| 14 | uxy_b             | xc, yc           | m/yr        | Basal horizontal velocity magnitude           | 1e-3,1e-4r |
| 15 | ux_s              | xc, yc           | m/yr        | Surface x-velocity                            | 1e-3,1e-4r |
| 16 | uy_s              | xc, yc           | m/yr        | Surface y-velocity                            | 1e-3,1e-4r |
| 17 | uz_s              | xc, yc           | m/yr        | Surface z-velocity                            | 1e-3,1e-4r |
| 18 | uxy_s             | xc, yc           | m/yr        | Surface horizontal velocity magnitude         | 1e-3,1e-4r |
| 19 | ux_i              | xc, yc, zeta     | m/yr        | Shearing x-velocity                           | 1e-3,1e-4r |
| 20 | uy_i              | xc, yc, zeta     | m/yr        | Shearing y-velocity                           | 1e-3,1e-4r |
| 21 | ux_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing x-velocity            | 1e-3,1e-4r |
| 22 | uy_i_bar          | xc, yc           | m/yr        | Depth-averaged shearing y-velocity            | 1e-3,1e-4r |
| 23 | uxy_i_bar         | xc, yc           | m/yr        | Depth-averaged horizontal velocity magnitude  | 1e-3,1e-4r |
| 24 | duxydt            | xc, yc           | m/yr^2      | Time derivative of uxy                        | 1e-4r      |
| 25 | duxdz             | xc, yc, zeta     | 1/yr        | x-velocity vertical gradient                  | 1e-4r      |
| 26 | duydz             | xc, yc, zeta     | 1/yr        | y-velocity vertical gradient                  | 1e-4r      |
| 27 | taud_acx          | xc, yc           | Pa          | Driving stress (x-dir)                        | 1e-4r      |
| 28 | taud_acy          | xc, yc           | Pa          | Driving stress (y-dir)                        | 1e-4r      |
| 29 | taud              | xc, yc           | Pa          | Driving stress magnitude                      | 1e-4r      |
| 30 | taub_acx          | xc, yc           | Pa          | Basal stress (x-dir)                          | 1e-4r      |
| 31 | taub_acy          | xc, yc           | Pa          | Basal stress (y-dir)                          | 1e-4r      |
| 32 | taub              | xc, yc           | Pa          | Basal stress magnitude                        | 1e-4r      |
| 33 | taul_int_acx      | xc, yc           | Pa          | Depth-integrated lateral stress (x-dir)       | 1e-4r      |
| 34 | taul_int_acy      | xc, yc           | Pa          | Depth-integrated lateral stress (y-dir)       | 1e-4r      |
| 35 | qq_gl_acx         | xc, yc           | m^3/yr      | Flux across grounding line                    | 1e-4r      |
| 36 | qq_gl_acy         | xc, yc           | m^3/yr      | Flux across grounding line                    | 1e-4r      |
| 37 | qq_acx            | xc, yc           | m^3/yr      | Flux (x-dir)                                  | 1e-4r      |
| 38 | qq_acy            | xc, yc           | m^3/yr      | Flux (y-dir)                                  | 1e-4r      |
| 39 | qq                | xc, yc           | m^3/yr      | Flux magnitude                                | 1e-4r      |
| 40 | visc_eff          | xc, yc, zeta     | Pa yr       | Effective viscosity                           | 1e-4r      |
| 41 | visc_eff_int      | xc, yc           | Pa yr m     | Depth-integrated viscosity                    | 1e-4r      |
| 42 | N_eff             | xc, yc           | Pa          | Effective pressure                            | 1e-4r      |
| 43 | cb_tgt            | xc, yc           | Pa          | Target basal parameter                        | 1e-4r      |
| 44 | cb_ref            | xc, yc           | --          | Reference basal parameter                     | 1e-4r      |
| 45 | c_bed             | xc, yc           | Pa          | Basal drag coefficient                        | 1e-4r      |
| 46 | f_slide           | xc, yc           | -           | Frozen-bed sliding-speed factor               | 1e-4       |
| 47 | beta_acx          | xc, yc           | Pa yr m^-1  | Basal stress factor (x)                       | 1e-4r      |
| 48 | beta_acy          | xc, yc           | Pa yr m^-1  | Basal stress factor (y)                       | 1e-4r      |
| 49 | beta              | xc, yc           | Pa yr m^-1  | Basal stress factor mag.                      | 1e-4r      |
| 50 | beta_eff          | xc, yc           | Pa yr m^-1  | Effective basal factor                        | 1e-4r      |
| 51 | f_vbvs            | xc, yc           | -           | Vertical basal stress                         | 1e-4       |
| 52 | ssa_mask_acx      | xc, yc           | -           | SSA mask (x-dir)                              | -          |
| 53 | ssa_mask_acy      | xc, yc           | -           | SSA mask (y-dir)                              | -          |
| 54 | ssa_err_acx       | xc, yc           | 1           | SSA error (x-dir)                             | 1e-4r      |
| 55 | ssa_err_acy       | xc, yc           | 1           | SSA error (y-dir)                             | 1e-4r      |
| 56 | jvel_dxx          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdx             | 1e-4r      |
| 57 | jvel_dxy          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdy             | 1e-4r      |
| 58 | jvel_dxz          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duxdz             | 1e-4r      |
| 59 | jvel_dyx          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydx             | 1e-4r      |
| 60 | jvel_dyy          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydy             | 1e-4r      |
| 61 | jvel_dyz          | xc, yc, zeta     | 1/yr        | Velocity Jacobian component duydz             | 1e-4r      |
| 62 | jvel_dzx          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdx             | 1e-4r      |
| 63 | jvel_dzy          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdy             | 1e-4r      |
| 64 | jvel_dzz          | xc, yc, zeta_ac  | 1/yr        | Velocity Jacobian component duzdz             | 1e-4r      |
| 65 | H_ice_solv        | xc, yc           | m           | Ice thickness of the last velocity solution   | 1e-3       |
| 66 | f_ice_solv        | xc, yc           | 1           | Ice fraction of the last velocity solution    | 1e-4       |
| 67 | uz_srf_err        | xc, yc           | m/yr        | Surface uz_star + smb (kinematic mismatch)    | 1e-4       |
