#!/bin/bash
#
# ISMIP7 Antarctica optimization spin-up (climate = esm). opt.cf_init=-1 starts
# cb_ref from the till friction of the bed (cb_tgt, from the ytill parameters).

resolution=ANT-16KM
output_path=output/ismip7/${resolution}/opt-l21-bedmap3

ctrl_params=(
    "domain.grid_hub=${resolution}"
    "isos.rheology_file=isostasy_data/earth_structure/yelmo/${resolution}_GIA_HR24.nc"
    "ctrl.run_step=spinup"
    "comps.surface_method=climate"
    "sim.relax=True"
    "sim.opt=True"
    "spinup.time_end=15.0e3"
    "sim.init_kill_shelves=True"
    "yelmo.nz_aa=11"
    "yelmo.dt_min=0.1"
    "tm_1D.dt=1.0"
    "tm_2Dsm.dt=2.5e3"
    "tm_2D.dt=15e3"
    "output.write_cmip=False"
)

opt_params=(
    "opt.H0=100"
    "opt.cf_time_end=15e3"
    "opt.tf_time_end=15e3"
    "opt.tau_c=500.0"
    "relax.tau1=100.0"
    "relax.time1=100.0"
    "relax.tau2=100.0"
    "relax.time2=100.0"
    "opt.cf_init=-1"
    "opt.H_grnd_lim=500.0"
    "ytill.scale_zb=1"
    "ytill.z0=-1000,-750,-500"
    "ytill.z1=0,500,1000"
    "ytill.cf_min=1e-3"
    "ytill.cf_ref=1e-0"
    "marine_shelf.gamma_quad_nl=14.5e3"
)

topo_params=(
    "ytopo.bmb_gl_method=pmp"
    "ytopo.gl_sep=3"
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

mat_params=(
    "ymat.enh_shear=1.0"
    "ymat.enh_stream=1.0"
    "ymat.enh_shlf=0.5"
)

runme -rs -q shared -w 2-00:00:00 -m 20G -e yelmox --omp 8 -n yelmox/yelmox_esm_Antarctica_ismip7.nml -o "${output_path}" \
      -p "${ctrl_params[@]}" "${opt_params[@]}" "${topo_params[@]}" "${calv_params[@]}" "${dyn_params[@]}" "${mat_params[@]}"
