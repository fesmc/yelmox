# ytopo

| id | variable          | dimensions  | units       | long_name                                          | prec  |
|----|-------------------|-------------|-------------|----------------------------------------------------|-------|
|  1 | H_ice             | xc, yc      | m           | Ice thickness                                      | 1e-3  |
|  2 | dHidt             | xc, yc      | m/yr        | Ice thickness rate of change                       | 1e-4  |
|  3 | dHidt_dyn         | xc, yc      | m/yr        | Ice thickness change due to dynamics               | 1e-4  |
|  4 | mb_net            | xc, yc      | m/yr        | Actual mass balance applied                        | 1e-4  |
|  5 | mb_relax          | xc, yc      | m/yr        | Change in mass balance due to relaxation           | 1e-4  |
|  6 | mb_resid          | xc, yc      | m/yr        | Residual mass balance                              | 1e-4  |
|  7 | mb_err            | xc, yc      | m/yr        | Residual error in mass balance accounting          | 1e-4  |
|  8 | smb               | xc, yc      | m/yr        | Surface mass balance                               | 1e-4  |
|  9 | bmb               | xc, yc      | m/yr        | Combined basal mass balance                        | 1e-4  |
| 10 | fmb               | xc, yc      | m/yr        | Combined frontal mass balance                      | 1e-4  |
| 11 | dmb               | xc, yc      | m/yr        | Subgrid discharge mass balance                     | 1e-4  |
| 12 | cmb               | xc, yc      | m/yr        | Calving mass balance                               | 1e-4  |
| 13 | bmb_ref           | xc, yc      | m/yr        | Reference basal mass balance                       | 1e-4  |
| 14 | fmb_ref           | xc, yc      | m/yr        | Reference frontal mass balance                     | 1e-4  |
| 15 | dmb_ref           | xc, yc      | m/yr        | Reference subgrid discharge mass balance           | 1e-4  |
| 16 | cmb_flt           | xc, yc      | m/yr        | Calving mass balance, floating cells               | 1e-4  |
| 17 | cmb_grnd          | xc, yc      | m/yr        | Calving mass balance, grounded cells               | 1e-4  |
| 18 | z_srf             | xc, yc      | m           | Surface elevation                                  | 1e-3  |
| 19 | dzsdt             | xc, yc      | m/yr        | Surface elevation rate of change                   | 1e-4  |
| 20 | eps_eff           | xc, yc      | 1/yr        | Effective strain rate                              | 1e-4r |
| 21 | tau_eff           | xc, yc      | Pa          | Effective stress                                   | 1e-4r |
| 22 | z_base            | xc, yc      | m           | Ice-base elevation                                 | 1e-3  |
| 23 | dzsdx             | xc, yc      | m/m         | Surface elevation slope, acx nodes                 | 1e-4r |
| 24 | dzsdy             | xc, yc      | m/m         | Surface elevation slope, acy nodes                 | 1e-4r |
| 25 | dHidx             | xc, yc      | m/m         | Ice thickness gradient, acx nodes                  | 1e-4r |
| 26 | dHidy             | xc, yc      | m/m         | Ice thickness gradient, acy nodes                  | 1e-4r |
| 27 | dzbdx             | xc, yc      | m/m         | Bedrock slope, acx nodes                           | 1e-4r |
| 28 | dzbdy             | xc, yc      | m/m         | Bedrock slope, acy nodes                           | 1e-4r |
| 29 | H_eff             | xc, yc      | m           | Effective ice thickness (margin-corrected)         | 1e-3  |
| 30 | H_grnd            | xc, yc      | m           | Grounded ice thickness                             | 1e-3  |
| 31 | H_calv            | xc, yc      | m           | Calving parameter field, ice thickness limit       | 1e-3  |
| 32 | kt_calv           | xc, yc      |             | Calving parameter field, vm-l19                    | 1e-4r |
| 33 | z_bed_filt        | xc, yc      | m           | Filtered bedrock elevation                         | 1e-3  |
| 34 | f_grnd            | xc, yc      |             | Grounded fraction                                  | 1e-4  |
| 35 | f_grnd_acx        | xc, yc      |             | Grounded fraction (acx nodes)                      | 1e-4  |
| 36 | f_grnd_acy        | xc, yc      |             | Grounded fraction (acy nodes)                      | 1e-4  |
| 37 | f_grnd_ab         | xc, yc      |             | Grounded fraction (ab nodes)                       | 1e-4  |
| 38 | f_ice             | xc, yc      |             | Ice-covered fraction                               | 1e-4  |
| 39 | f_grnd_bmb        | xc, yc      |             | Grounded fraction for basal mass balance           | 1e-4  |
| 40 | f_grnd_pin        | xc, yc      |             | Grounded fraction from subgrid pinning points      | 1e-4  |
| 41 | dist_margin       | xc, yc      | m           | Distance to nearest margin point                   | 1e-3  |
| 42 | dist_grline       | xc, yc      | m           | Distance to nearest grounding-line point           | 1e-3  |
| 43 | mask_bed          | xc, yc      |             | Multi-valued bed mask                              | -     |
| 44 | mask_grz          | xc, yc      |             | Multi-valued grounding-line zone mask              | -     |
| 45 | mask_frnt         | xc, yc      |             | Multi-valued ice front mask                        | -     |
| 46 | dHidt_dyn_raw_n   | xc, yc      | m/yr        | Ice thickness change due to advection (previous)   | 1e-4  |
| 47 | H_ice_n           | xc, yc      | m           | Ice thickness from previous timestep               | 1e-3  |
| 48 | z_srf_n           | xc, yc      | m           | Surface elevation from previous timestep           | 1e-3  |
| 49 | H_ice_dyn         | xc, yc      | m           | Dynamic ice thickness                              | 1e-3  |
| 50 | f_ice_dyn         | xc, yc      |             | Dynamic ice-covered fraction                       | 1e-4  |
| 51 | pc_pred_H_ice     | xc, yc      | m           | Predicted ice thickness                            | 1e-3  |
| 52 | pc_pred_dHidt_dyn | xc, yc      | m/yr        | Predicted dynamic ice thickness rate of change     | 1e-4  |
| 53 | pc_pred_mb_net    | xc, yc      | m/yr        | Predicted net mass balance                         | 1e-4  |
| 54 | pc_pred_smb       | xc, yc      | m/yr        | Predicted surface mass balance                     | 1e-4  |
| 55 | pc_pred_bmb       | xc, yc      | m/yr        | Predicted basal mass balance                       | 1e-4  |
| 56 | pc_pred_fmb       | xc, yc      | m/yr        | Predicted frontal mass balance                     | 1e-4  |
| 57 | pc_pred_dmb       | xc, yc      | m/yr        | Predicted discharge mass balance                   | 1e-4  |
| 58 | pc_pred_cmb       | xc, yc      | m/yr        | Predicted calving mass balance                     | 1e-4  |
| 59 | pc_corr_H_ice     | xc, yc      | m           | Corrected ice thickness                            | 1e-3  |
| 60 | pc_corr_dHidt_dyn | xc, yc      | m/yr        | Corrected dynamic ice thickness rate of change     | 1e-4  |
| 61 | pc_corr_mb_net    | xc, yc      | m/yr        | Corrected net mass balance                         | 1e-4  |
| 62 | pc_corr_smb       | xc, yc      | m/yr        | Corrected surface mass balance                     | 1e-4  |
| 63 | pc_corr_bmb       | xc, yc      | m/yr        | Corrected basal mass balance                       | 1e-4  |
| 64 | pc_corr_fmb       | xc, yc      | m/yr        | Corrected frontal mass balance                     | 1e-4  |
| 65 | pc_corr_dmb       | xc, yc      | m/yr        | Corrected discharge mass balance                   | 1e-4  |
| 66 | pc_corr_cmb       | xc, yc      | m/yr        | Corrected calving mass balance                     | 1e-4  |
| 67 | tau_relax         | xc, yc      |             | Relaxed tau time                                   | 1e-4r |
| 68 | lsf               | xc, yc      |             | Level-set function                                 | 1e-4r |
| 69 | cmb_flt_x         | xc, yc      | m/yr        | Floating calving rate (x-direction)                | 1e-4  |
| 70 | cmb_flt_y         | xc, yc      | m/yr        | Floating calving rate (y-direction)                | 1e-4  |
| 71 | dzsdt_kin         | xc, yc      | m/yr        | Kinematic rate of the ice-column surface           | 1e-4  |
| 72 | dzbdt_kin         | xc, yc      | m/yr        | Kinematic rate of the ice-column base              | 1e-4  |
| 73 | calv_rate_flt     | xc, yc      | m/yr        | Calving speed of floating front cells              | 1e-4  |
| 74 | calv_rate_grnd    | xc, yc      | m/yr        | Calving speed of grounded front cells              | 1e-4  |
| 75 | mb_clip           | xc, yc      | m/yr        | Clip of negative ice thickness after transport     | 1e-4  |
| 76 | mask_kin          | xc, yc      |             | Column rate given by the applied thickness step    | -     |
| 77 | mb_border         | xc, yc      | m/yr        | Reset of outflow-border cells (domain outflow)     | 1e-4  |
| 78 | tau_vm            | xc, yc      | Pa          | Tensile von Mises stress (Morlighem et al., 2016)  | 1e-4r |
