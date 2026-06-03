#!/bin/bash

resolution=32km
output_path=output_albedo/ismip7/opt-${resolution}_energy_cbtgtVA26_enh3_jbn_smb

ctrl_params=(
    "ctrl.run_step=spinup"
    "esm.use_smb=True"
    "spinup.equil_method=opt"
    "spinup.time_end=15.0e3"
    "spinup.kill_shelves=True"
    "tm_1D.dt=10.0"
    "tm_2Dsm.dt=2e3"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
)

opt_params=(
    "opt.H0=-1"
	"opt.cf_time_end=15e3"
    "opt.tf_time_end=15e3"
    "opt.tau_c=100.0"
    "opt.rel_tau1=1.0"
    "opt.rel_time1=10.0"
    "opt.rel_tau2=1.0"
    "opt.rel_time2=10.0"
	"opt.opt_cf_min=1e-3"
    "opt.cf_init=-1"
	"opt.basin_fill=True"
    "ytill.scale_zb=1"
    "ytill.z0=-700"
    "ytill.z1=700"
    "ytill.cf_min=0.1"
    "ytill.cf_ref=0.4"
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
    "ymat.enh_shear=3.0"
    "ymat.enh_stream=3.0"
    "ymat.enh_shlf=1.0"
)

runme -rs -q 48h -e esm -n par/yelmo_Antarctica_esm_ismip7.nml -o "${output_path}" \
      -p "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${neff_params[@]}" "${mat_params[@]}"
