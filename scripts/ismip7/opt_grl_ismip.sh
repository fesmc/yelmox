#!/bin/bash
#
# ISMIP7 Greenland optimization spin-up (L. Gutierrez Gonzalez), ported from
# legacy/opt_grl_ismip.sh (yelmox v2.2, par/yelmo_Greenland_esm_ismip7.nml) to
# yelmox (climate = esm) and yelmox/yelmox_esm_Greenland.nml. Only the
# overrides of the original are carried over; the base configuration is the
# current par file. Key changes:
#   spinup.equil_method       -> coupling.equil_method
#   opt.opt_cf_min=0.002      -> ytill.cf_min=0.002 (one lower bound for the
#                                optimization and cb_tgt; was ytill.cf_min=1e-1)
#   opt.use_yelmo_cf_min      -> removed
#   fhyd.*                    -> yhyd.* (the original set bkt_N_closure and
#                                till_delta twice; runme keeps the last values)
#   ycalv.tau_ice             -> ycalv.tau_ice_flt
#   ydyn.scale_T=0            -> ydyn.slide_T=False
#   gcm_to_ref/gcm_so_ref     -> dropped: input/esm/esm_grl_ismip7.nml already
#                                uses the annual 1981-2010 mean (data_params_8KM)
# opt.cf_init=-1 starts cb_ref from the till friction of the bed (cb_tgt).

resolution=GRL-8KM
output_path=output/ismip7/${resolution}/spinup/26_tf1

res_params=(
    "domain.grid_hub=${resolution}"
)

ctrl_params=(
    "ctrl.run_step=spinup"
    "coupling.smb_method=climate"
    "coupling.equil_method=opt"
    "spinup.time_end=15e3"
    "coupling.kill_shelves=True"
    "tm_1D.dt=1.0"
    "tm_2Dsm.dt=500"
    "tm_2D.dt=15e3"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
    "output.write_cmip=False"
)

opt_params=(
    "opt.H0=100"
    "opt.sigma_err=50e3"
    "opt.cf_time_end=15e3"
    "opt.tf_time_end=15e3"
    "opt.tau_c=500"
    "opt.rel_tau1=100.0"
    "opt.rel_time1=100.0"
    "opt.rel_tau2=100.0"
    "opt.rel_time2=100.0"
    "opt.sigma_vel=100"
    "opt.cf_init=-1"
    "opt.H_grnd_lim=500.0"
    "ytill.scale_zb=1"
    "ytill.z0=-500"
    "ytill.z1=500"
    "ytill.cf_min=0.002"
    "ytill.cf_ref=1e-0"
)

hyd_params=(
    "yhyd.bkt_N_closure=2"
    "yhyd.till_delta=0.4"
)

topo_params=(
    "ytopo.bmb_gl_method=pmp"
    "ytopo.gl_sep=2"
    "ytopo.fmb_method=1"
    "ytopo.fmb_scale=10"
)

calv_params=(
    "ycalv.calv_flt_method=vm-l19"
    "ycalv.calv_grnd_method=zero"
    "ycalv.tau_ice_flt=200e3"
    "ycalv.kt_ref=0.003"
    "ycalv.use_lsf=False"
    "marine_shelf.bmb_method=lin"
    "marine_shelf.tf_method=1"
    "marine_shelf.gamma_lin=300"
)

dyn_params=(
    "ydyn.beta_min=50"
    "ydyn.solver=diva"
    "ydyn.slide_T=False"
    "ydyn.ssa_solver=energy"
    "ydyn.ssa_lat_bc=all"
)

mat_params=(
    "ymat.enh_shear=1.0"
    "ymat.enh_stream=1.0"
    "ymat.enh_shlf=0.5"
)

runme -rs -q 48h -e yelmox --omp 8 -n yelmox/yelmox_esm_Greenland.nml -o "${output_path}" \
      -p "${res_params[@]}" "${ctrl_params[@]}" "${opt_params[@]}" "${hyd_params[@]}" \
      "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${mat_params[@]}"
