# ybound

| id | variable          | dimensions  | units       | long_name                                              | prec  |
|----|-------------------|-------------|-------------|--------------------------------------------------------|-------|
|  1 | z_bed             | xc, yc      | m           | Bedrock elevation                                      | 1e-3  |
|  2 | z_bed_sd          | xc, yc      | m           | Standard deviation of bedrock elevation                | 1e-3  |
|  3 | z_sl              | xc, yc      | m           | Sea level elevation                                    | 1e-3  |
|  4 | H_sed             | xc, yc      | m           | Sediment thickness                                     | 1e-3  |
|  5 | smb_ref           | xc, yc      | m/yr        | Surface mass balance                                   | 1e-4  |
|  6 | T_srf             | xc, yc      | K           | Surface temperature                                    | 1e-3  |
|  7 | bmb_shlf          | xc, yc      | m/yr        | Basal mass balance for ice shelf                       | 1e-4  |
|  8 | fmb_shlf          | xc, yc      | m/yr        | Frontal mass balance for ice shelf                     | 1e-4  |
|  9 | T_shlf            | xc, yc      | K           | Ice shelf temperature                                  | 1e-3  |
| 10 | Q_geo             | xc, yc      | mW m^-2     | Geothermal heat flow at depth                          | 1e-4r |
| 11 | enh_srf           | xc, yc      | -           | Enhancement factor at the surface                      | 1e-4r |
| 12 | basins            | xc, yc      | -           | Basin ids (0 = none)                                   | -     |
| 13 | basin_mask        | xc, yc      | -           | Mask for basins                                        | -     |
| 14 | regions           | xc, yc      | -           | Region codes, deepest level (0 = none)                 | -     |
| 15 | calv_mask         | xc, yc      | -           | Locations where calving is not allowed                 | -     |
| 16 | H_ice_ref         | xc, yc      | m           | Reference ice thickness for relaxation routines        | 1e-3  |
| 17 | z_bed_ref         | xc, yc      | m           | Reference bedrock elevation for relaxation routines    | 1e-3  |
| 18 | mask_ice          | xc, yc      | -           | Ice mask (0=none, 1=fixed, 2=dynamic)                  | -     |
| 19 | z_bed_n           | xc, yc      | m           | Bedrock elevation at the previous yelmo_update call    | 1e-3  |
| 20 | z_sl_n            | xc, yc      | m           | Sea level at the previous yelmo_update call            | 1e-3  |
