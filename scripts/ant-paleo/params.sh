# Antarctic paleo setup on the base yelmox/yelmox_Antarctica.nml, as runme -p
# arrays. Sourced by the run_*.sh scripts of this folder and by
# scripts/docs-figures/run.sh. The named groups (series_*, term_*, snapshot_*)
# come from input/forcing_db.nml (docs/_forcing-terms.md).

# 32 km physics and data of the spin-up and the transient.
antpal_physics=(
    "domain.grid_hub=ANT-32KM"          # single-grid setup: the hub on the Yelmo grid
    "domain.grid_shelf=ANT-32KM"        # marine shelf on the Yelmo grid
    "ctrl.dtt=5.0"                      # [yr] main loop timestep
    "forcing.series==dTo"               # dTo: the mean sea-level air anomaly x series_dTo.factor (derived)
    "forcing.terms==ocean"              # the ocean profiles plus dTo
    # surface: chion with the ITM (group [itm], shared with smbpal)
    "chion.model=itm"                   # snowpack model (parameters: [itm], the chion default nml_itm)
    "smbpal.abl_method=itm"             # smbpal ablation
    "itm.itm_c=-55.0"                   # ITM melt parameters
    "itm.itm_b=3.0"
    "itm.itm_lat0=-60.0"
    "itm.alb_ice=0.7"                   # ice albedo
    "yelmo_data.pd_tsrf_path=ice_data/{domain}/{grid_name}/{grid_name}_RACMO23-ERAINT-HYBRID_1981-2010.nc"
    "yelmo_data.pd_smb_path=ice_data/{domain}/{grid_name}/{grid_name}_RACMO23-ERA-INTERIM_monthly_1981-2010.nc"
    # solid Earth and sea level
    "isos.pad=512."                     # [km] padding around the domain
    "isos.viscosity_scaling_method=stddev"
    "isos.viscosity_scaling=0"
    "isos.tau=2000.0"                   # [yr] relaxation time (method 2)
    "barysealevel.sl_path=input/sealevel_rohling_450kyr.dat"   # sea level series (method = file)
    "ghf.obs_path=fesmdata/{domain}/{grid_name}/{grid_name}_GHF-Loesing2021.nc"   # geothermal heat flux
    "ghf.obs_err_name=None"             # ... without an error field
    # ice dynamics
    "ytopo.grad_lim=0.1"                # [m/m] maximum slope in gradients
    "ytopo.gz_Hg1=50"                   # grounding zone: penetration limit of bmb_shlf
    "ytopo.fmb_method=1"                # front melt proportional to bmb_shlf
    "ycalv.calv_flt_method=vm-m16"      # von Mises calving, floating fronts
    "ycalv.calv_grnd_method=zero"       # no calving of grounded fronts
    "ycalv.k2=3e9"                      # [m yr] eigen calving factor
    "ycalv.w2=25"                       # weight of the 2nd principal stress
    "ycalv.tau_ice_flt=1500.0e3"        # [Pa] ice strength, floating fronts (vm-m16)
    "ytill.method=-1"                   # cb_ref set externally (optimization)
    "ytill.scale_zb=2"                  # exponential scaling with bed elevation
    "ytill.z0=-200.0"                   # [m] its lower limit
    "ytill.cf_min=0.001"                # minimum cf
    "ytill.cf_ref=0.05"                 # reference cf
    "yhyd.bkt_N_closure=3"              # effective pressure: till
    "yhyd.till_delta=0.8"
    "ymat.enh_shlf=0.7"                 # enhancement factor, shelves
    # thermal forcing corrections by Zwally basin ([marine_shelf] corr_b<k>)
    "marine_shelf.corr_names==b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18"
    "marine_shelf.corr_values==0.0,0.0,0.0,0.0,0.0,0.9,0.0,0.0,0.0,0.0,0.0,0.5,0.0,-1.0,-1.0,-0.7,-0.7,0.0"
    # optimization and relaxation of the spin-up
    "opt.cf_time_end=7.5e3"             # [yr] end of the cf_ref optimization
    "opt.cf_init=-1.0"                  # initial cf_ref from cb_tgt
    "opt.fill_method=target"            # cf_ref outside the ice
    "opt.tf_time_end=7.5e3"             # [yr] end of the tf_corr optimization
    "relax.time1=1000"                  # [yr] end of the tau1 period
    "relax.time2=3000"                  # [yr] end of the relaxation
)

# Last glacial period forcing (run_lgp.sh): PMIP3 piControl and LGM blended by
# the glacial index input/alpha_combined_125kyr_interp.dat; Antarctic summer
# (DJF) as the second temperature of the PMIP3 snapshots.
antpal_lgp=(
    "forcing.series==alpha,dTo"         # alpha: the glacial index ([series_alpha] of the base par file); dTo as in antpal_physics
    "forcing.terms==lgm,ocean"          # lgm: lgm - picontrol (PMIP3) weighted by alpha, on tas and pr
    "snapshot_picontrol.tas==pmip3_mean_tas_ann,pmip3_mean_tas_djf"   # summer: DJF
    "snapshot_lgm.tas==pmip3_mean_tas_ann,pmip3_mean_tas_djf"         # summer: DJF
    "orbit.method=forcing_time"                # transient insolation (orbit at the forcing time)
    "barysealevel.method=file"          # sea level from sl_path
)
