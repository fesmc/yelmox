#!/bin/bash

resolution=GRL-8KM
temp_path=/p/projects/megarun/luciagu/yelmox_v2.2/ismip7_output/${resolution}/spinup/26_tf1

res_params=(
    "yelmo.grid_name=${resolution}"
)

ctrl_params=(
    "ctrl.run_step=spinup"
    "esm.use_smb=True"
    "spinup.equil_method=opt"
    "spinup.time_end=15e3"
    "spinup.kill_shelves=True"
    "tm_1D.dt=1.0"
    "tm_2Dsm.dt=500"
    "tm_2D.dt=15e3"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
    "esm.write_formatted=False"
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
    "opt.use_yelmo_cf_min=False"
    "opt.opt_cf_min=0.002"
    "opt.sigma_vel=100"
    "opt.cf_init=-1"
    "opt.H_grnd_lim=500.0"
    "ytill.scale_zb=1"
    "ytill.z0=-500"
    "ytill.z1=500"
    "ytill.cf_min=1e-1"
    "ytill.cf_ref=1e-0"
    "fhyd.bkt_N_closure=3"
    "fhyd.till_delta=0.04"
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
    "ycalv.tau_ice=200e3"
    "marine_shelf.bmb_method=lin"
    "marine_shelf.tf_method=1"
    "marine_shelf.gamma_lin=300"
    "ycalv.kt_ref=0.003"
    "ycalv.use_lsf=False"
    "fhyd.till_delta=0.4"
    "fhyd.bkt_N_closure=2"
)

dyn_params=(
    "ydyn.beta_min=50"
    "ydyn.solver=diva"
    "ydyn.scale_T=0"
	"ydyn.ssa_solver=energy"
	"ydyn.ssa_lat_bc=all"
)
  

mat_params=(
    "ymat.enh_shear=1.0"
    "ymat.enh_stream=1.0"
    "ymat.enh_shlf=0.5"
)

data_params_8KM=(
    "gcm_to_ref.with_time=False"
    "gcm_to_ref.time_par=1981 2010 0 1"
    "gcm_so_ref.with_time=False"
    "gcm_so_ref.time_par=1981 2010 0 1"
)

data_params_16KM=(
    "gcm_to_ref.with_time=True"
    "gcm_to_ref.time_par=1981 2010 0 12"
    "gcm_so_ref.with_time=True"
    "gcm_so_ref.time_par=1981 2010 0 12"
)

nohup runme -rs -q standby -w 24:00:00 -e esm --omp 32 -n par/yelmo_Greenland_esm_ismip7.nml -o "${temp_path}" \
      -p "${res_params[@]}" "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" \
      "${dyn_params[@]}" "${mat_params[@]}" "${data_params_8KM[@]}"

