#!/bin/bash

resolution=16km
output_path=output_albedo/1pctCO2/opt-${resolution}-l21-bedmap3

ctrl_params=(
    "ctrl.run_step=spinup"
    "esm.use_smb=False"
    "spinup.equil_method=opt"
    "spinup.time_end=15.0e3"
    "spinup.kill_shelves=True"
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
    "opt.use_yelmo_cf_min=False"
    "opt.opt_cf_min=1e-3"
    "opt.cf_init=-1"
    "opt.H_grnd_lim=500.0"
    "ytill.scale_zb=1"
    "ytill.z0=-1000,-500"
    "ytill.z1=0,500"
    "ytill.cf_min=1e-1"
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
    "ycalv.tau_ice=200e3"
)

dyn_params=(
    "ydyn.beta_min=10.0"
    "ydyn.solver=diva"
    "ydyn.scale_T=0"
    "ydyn.ssa_solver=energy"
    "ydyn.ssa_lat_bc=all"
)

neff_params=(
    "yneff.method=2"
    "yneff.p=1.0"
)     

mat_params=(
    "ymat.enh_shear=1.0"
    "ymat.enh_stream=1.0"
    "ymat.enh_shlf=0.5"
)

runme -rs -q 48h -e esm --omp 8 -n par/yelmo_Antarctica_esm_ismip7.nml -o "${output_path}" \
      -p "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${neff_params[@]}" "${mat_params[@]}"
