#!/bin/bash

resolution=ANT-16KM
output_path=output_albedo/1pctCO2/opt-${resolution}-l21-bedmap3

ctrl_params=(
    "domain.grid_hub=${resolution}"
    "ctrl.run_step=spinup"
    "comps.smb_method=smbpal"
    "sim.opt=True"
    "spinup.time_end=15.0e3"
    "sim.init_kill_shelves=True"
    "tm_1D.dt=10.0"
    "tm_2Dsm.dt=2e3"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
    "marine_shelf.extrap_shlf=True"
)

opt_params=(
    "opt.H0=100"
    "opt.cf_time_end=20e3"
    "opt.tf_time_end=20e3"
    "opt.tau_c=500.0"
    "opt.rel_tau1=100.0"
    "opt.rel_time1=100.0"
    "opt.rel_tau2=100.0"
    "opt.rel_time2=100.0"
    "opt.cf_init=-1"
    "opt.H_grnd_lim=500.0"
    "ytill.scale_zb=1"
    "ytill.z0=-1000,-500"
    "ytill.z1=0,500"
    "ytill.cf_min=1e-3"
    "ytill.cf_ref=1e-0"
    "marine_shelf.gamma_quad_nl=14.5e3"
)

topo_params=(
    "ytopo.bmb_gl_method=pmp"
    "ytopo.gl_sep=2"
)

calv_params=(
    "ycalv.calv_flt_method=equil"
    "ycalv.calv_grnd_method=equil"
    "ycalv.tau_ice_flt=200e3"
)

dyn_params=(
    "ydyn.beta_min=10.0"
    "ydyn.solver=diva"
    "ydyn.ssa_solver=energy"
    "ydyn.ssa_lat_bc=all"
)

hyd_params=(
    "yhyd.bkt_N_closure=2"
    "yhyd.marine_p=1.0"
)     

mat_params=(
    "ymat.enh_shear=1.0"
    "ymat.enh_stream=1.0"
    "ymat.enh_shlf=0.5"
)

runme -rs -q 48h -e yelmox --omp 8 -n yelmox/yelmox_esm_Antarctica_1pctCO2.nml -o "${output_path}" \
      -p "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${hyd_params[@]}" "${mat_params[@]}"
