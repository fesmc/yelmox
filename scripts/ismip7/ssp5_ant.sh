#!/bin/bash

resolution=32km
output_path=output_albedo/ismip7/ssp585-${resolution}_test_formatted


ctrl_params=(
    "ctrl.run_step=transient"
    "esm.use_esm=True"
    "esm.use_smb=True"
    "esm.use_proj=True"
    "transient.equil_method=none"
    "transient.time_init=2016"
    "transient.time_end=2300"
    "spinup.kill_shelves=False"
    "tm_1D.dt=1.0"
    "tm_2Dsm.method=const"
    "tm_2Dsm.dt=10.0"
    "tm_2D.method=times"
    "tm_2D.times=2016,2300"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
    "esm.write_formatted=True"
    "esm.dt_formatted=10"
)

topo_params=(
    "ytopo.bmb_gl_method=pmp"
    "ytopo.gl_sep=2"
)

calv_params=(
    "ycalv.calv_flt_method=vm-m16"
    "ycalv.calv_grnd_method=vm-m16"
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

path_restart=/albedo/home/jablas001/yelmo-awi/yelmox/output_albedo/ismip7/opt-${resolution}-l21-bedmap3_ensemble/3/restart-0.000-kyr/
#path_restart=/albedo/home/jablas001/yelmo-awi/yelmox/output_albedo/ismip7/opt-${resolution}-l21-bedmap3_ensemble/3/restart-0.000-kyr/
restart_params=(
  "yelmo.restart=${path_restart}/yelmo_restart.nc" 
  "marine_shelf.restart=${path_restart}/marine_shelf.nc" 
  "isos.restart=${path_restart}/isos_restart.nc" 
  "barysealevel.restart=${path_restart}/bsl_restart.nc"
)

runme -rs -q 48h -e esm --omp 8 -n par/yelmo_Antarctica_esm_ismip7.nml -o "${output_path}" \
      -p "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${neff_params[@]}" "${mat_params[@]}" "${restart_params[@]}"
