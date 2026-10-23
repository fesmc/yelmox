# ymat

| id | variable          | dimensions       | units       | long_name                                     | prec  |
|----|-------------------|------------------|-------------|-----------------------------------------------|-------|
|  1 | enh               | xc, yc, zeta     | -           | Enhancement factor                            | 1e-4r |
|  2 | enh_bnd           | xc, yc, zeta     | -           | Imposed enhancement factor                    | 1e-4r |
|  3 | enh_bar           | xc, yc           | -           | Depth-averaged enhancement                    | 1e-4r |
|  4 | ATT               | xc, yc, zeta     | -           | Rate factor                                   | 1e-4r |
|  5 | ATT_bar           | xc, yc           | -           | Depth-averaged rate factor                    | 1e-4r |
|  6 | visc              | xc, yc, zeta     | Pa yr       | Ice viscosity                                 | 1e-4r |
|  7 | visc_bar          | xc, yc           | Pa yr       | Depth-averaged ice viscosity                  | 1e-4r |
|  8 | visc_int          | xc, yc           | Pa yr m     | Ice viscosity interpolated at interfaces      | 1e-4r |
|  9 | f_shear_bar       | xc, yc           | -           | Depth-averaged shear fraction                 | 1e-4  |
| 10 | strn2D_dxx        | xc, yc           | 1/yr        | 2D strain rate tensor component dxx           | 1e-4r |
| 11 | strn2D_dyy        | xc, yc           | 1/yr        | 2D strain rate tensor component dyy           | 1e-4r |
| 12 | strn2D_dxy        | xc, yc           | 1/yr        | 2D strain rate tensor component dxy           | 1e-4r |
| 13 | strn2D_dxz        | xc, yc           | 1/yr        | 2D strain rate tensor component dxz           | 1e-4r |
| 14 | strn2D_dyz        | xc, yc           | 1/yr        | 2D strain rate tensor component dyz           | 1e-4r |
| 15 | strn2D_de         | xc, yc           | 1/yr        | 2D effective strain rate                      | 1e-4r |
| 16 | strn2D_div        | xc, yc           | 1/yr        | 2D horizontal divergence                      | 1e-4r |
| 17 | strn2D_f_shear    | xc, yc           |             | 2D strain rate shear fraction                 | 1e-4  |
| 18 | strn_dxx          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxx              | 1e-4r |
| 19 | strn_dyy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyy              | 1e-4r |
| 20 | strn_dxy          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxy              | 1e-4r |
| 21 | strn_dxz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dxz              | 1e-4r |
| 22 | strn_dyz          | xc, yc, zeta     | 1/yr        | Strain rate tensor component dyz              | 1e-4r |
| 23 | strn_de           | xc, yc, zeta     | 1/yr        | Effective strain rate                         | 1e-4r |
| 24 | strn_div          | xc, yc, zeta     | 1/yr        | Horizontal divergence                         | 1e-4r |
| 25 | strn_f_shear      | xc, yc, zeta     |             | Strain rate shear fraction                    | 1e-4  |
| 26 | strs2D_txx        | xc, yc           | Pa          | 2D stress tensor component txx                | 1e-4r |
| 27 | strs2D_tyy        | xc, yc           | Pa          | 2D stress tensor component tyy                | 1e-4r |
| 28 | strs2D_txy        | xc, yc           | Pa          | 2D stress tensor component txy                | 1e-4r |
| 29 | strs2D_txz        | xc, yc           | Pa          | 2D stress tensor component txz                | 1e-4r |
| 30 | strs2D_tyz        | xc, yc           | Pa          | 2D stress tensor component tyz                | 1e-4r |
| 31 | strs2D_te         | xc, yc           | Pa          | 2D effective stress                           | 1e-4r |
| 32 | strs2D_tau_eig_1  | xc, yc           | Pa          | 2D stress first principal eigenvalue          | 1e-4r |
| 33 | strs2D_tau_eig_2  | xc, yc           | Pa          | 2D stress second principal eigenvalue         | 1e-4r |
