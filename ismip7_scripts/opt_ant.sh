#!/bin/bash

resolution=32km
output_path=output_albedo/ismip7/opt-${resolution}_pmp_leguy

        main_params=(
            "ctrl.run_step=spinup"
            "spinup.equil_method=opt"
            "spinup.time_end=15.0e3"
	    "spinup.kill_shelves=True"
            "tm_1D.dt=10.0"
            "tm_2D.method=const"
            "tm_2Dsm.dt=2e3"
            "yelmo.nz_aa=11"
            "yelmo.dt_min=0.1"
            "ytopo.bmb_gl_method=pmp"
            "ytopo.gl_sep=2"
            "ydyn.beta_min=10.0"
            "ydyn.solver=diva"
            "yneff.method=2"
            "yneff.p=1.0"
	    "opt.H0=-1"
	    "opt.cf_time_end=15e3"
            "opt.tf_time_end=15e3"
	    "opt.tau_c=100.0"
            "opt.rel_tau1=1.0"
            "opt.rel_time1=10.0"
            "opt.rel_tau2=1.0"
            "opt.rel_time2=10.0"
	    "opt.cf_min=5e-4"
            "opt.cf_init=1e-2"
            "opt.tf_min=-1.0"
	    "opt.tf_max=1.0"
	    "opt.basin_fill=True"
            "ytill.cf_ref=1e-1"
	    "ycalv.calv_flt_method=equil"
            "ycalv.calv_grnd_method=equil"
            "ycalv.tau_ice=200e3"
	)

./runme -rs -q 48h -e esm -n par/ismip7/yelmo_Antarctica_esm_ismip7.nml -o "${output_path}" -p "${main_params[@]}"
