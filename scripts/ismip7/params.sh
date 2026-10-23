# ISMIP7 Antarctic setup on the base yelmox/yelmox_esm_Antarctica.nml, as a runme
# -p array. Sourced by run_ismip7_antarctica.sh and opt_ant.sh of this folder, by
# scripts/1pctCO2 (run_1pctco2_antarctica.sh, opt_ant.sh) and by
# scripts/docs-figures/run.sh.

# 16 km, the ISMIP7 observations (bedmap3 topography, Mouginot basins, geothermal
# heat flux, velocity), RACMO2.3p2 smb and T_srf, and the ISMIP7 ESM forcing
# (input/esm/esm_ant_ismip7.nml) with its smb.
ismip7_ant=(
    # time
    "spinup.dtt=1.0"                        # [yr] main loop timestep
    "spinup.time_hist==1850,2014"           # [yr] historical period
    "spinup.time_proj==2015,2300"           # [yr] projection period
    "spinup.time_ref==1985,2014"            # [yr] reference (RACMO2.3p2 climatology 1985-2014)
    "spinup.time_esm_ref==1960,1989"        # [yr] reference of the ESM
    "transient.tstep_const=2015.0"          # [yr] time of the "const" method
    "transient.time_init=2015.0"            # [yr CE]
    "transient.time_end=2300.0"             # [yr CE]
    "transient.time_hist==1850,2014"        # as [spinup]
    "transient.time_proj==2015,2300"
    "transient.time_ref==1985,2014"
    "transient.time_esm_ref==1960,1989"
    # components and grid
    "comps.surface_method=climate"          # the ESM smb
    "sim.relax=False"                       # no relaxation towards the reference topography
    "domain.grid_hub=ANT-16KM"              # hub, and through it every component grid
    # ISMIP7 observations: topography and basins on the hub grid (not FesmData)
    "topo.path=ice_data/ISMIP7/{domain}/{grid_name}/obs/{grid_name}_ObsISMIP7-v1.1.nc"
    "topo.vars==z_bed,H_ice,z_srf"
    "topo.names==bed_bedmap3,thickness_bedmap3,surface_bedmap3"   # variables of the file
    "regions.basin_sets="                   # no FesmData basin sets: the basins are the layer ismip7
    "regions.layers=ismip7"
    "regions.path_ismip7=ice_data/ISMIP7/{domain}/{grid_name}/obs/{grid_name}_ObsISMIP7-v1.1.nc"
    "regions.var_ismip7=mouginot_basins"
    "domain.basins=ismip7"                  # basins of the hub output
    "yelmo_masks.basins=ismip7"             # basins of Yelmo
    "esm.basins=ismip7"                     # basins of the ocean extrapolation
    "marine_shelf.basins=ismip7"            # basins of the marine shelf
    "marine_shelf.corr_b12==ismip7:12"      # basin of the correction b12 (== : the ":" is not a range)
    "ghf.obs_path=ice_data/ISMIP7/{domain}/{grid_name}/obs/{grid_name}_ObsISMIP7-v1.1.nc"
    "ghf.obs_name=geothermal_heat_flux2"    # geothermal heat flux
    "ghf.obs_err_name=none"                 # ... without an error field
    "yelmo_data.pd_vel_path=ice_data/ISMIP7/{domain}/{grid_name}/obs/{grid_name}_ObsISMIP7-v1.1.nc"
    "yelmo_data.pd_vel_names==vx_mosaic,vy_mosaic"   # surface velocity
    # present-day smb and T_srf: RACMO2.3p2 1981-2010
    "yelmo_data.pd_tsrf_path=fesmdata/{domain}/{grid_name}/{grid_name}_SMB-RACMO2.3p2-1981-2010.nc"
    "yelmo_data.pd_smb_path=fesmdata/{domain}/{grid_name}/{grid_name}_SMB-RACMO2.3p2-1981-2010.nc"
    # ESM forcing
    "esm.par_file=input/esm/esm_ant_ismip7.nml"
    "esm.esm_name=CESM2-WACCM"              # ESM of the forcing files
)
