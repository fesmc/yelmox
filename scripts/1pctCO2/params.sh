# 1pctCO2 setup as runme -p arrays. Sourced by the scripts of this folder and by
# scripts/docs-figures/run.sh. The forcing is the GCM's absolute tas/pr/thetao/so,
# self-referenced to the run start (~1xCO2): see the headers of
# input/esm/esm_grl_1pctCO2.nml and esm_ant_1pctCO2.nml. The ESM gives no smb, so
# the surface is chion.

# Greenland, on the base yelmox/yelmox_esm_Greenland.nml.
pct_grl=(
    "transient.tstep_const=2020"            # [yr] time of the "const" method
    "transient.time_init=2020.0"            # [yr] 1pctCO2 start (year 2020 = 1xCO2)
    "transient.time_end=2160.0"             # [yr] 1pctCO2 end (140 yr; ~4xCO2 near yr 140)
    "transient.time_ref==1961,1990"         # [yr] reference period of the MAR base climatology
    "transient.time_esm_ref==2020,2021"     # [yr] self-reference: dts = tas(t) - mean(tas over 2020-2021)
    "transient.time_hist==2018,2019"        # [yr] no historical period: kept below time_init (use_hist=False)
    "transient.time_proj==2020,2160"        # [yr] projection window (full run)
    "comps.surface_method=chion"            # the ESM gives no smb
    "esm.par_file=input/esm/esm_grl_1pctCO2.nml"
    "esm.experiment=1pctCO2"                # the spin-up sets ctrl
    "esm.esm_name=MPI-ESM1-2-LR"            # {gcm} of the forcing files
)

# Antarctica, on the base yelmox/yelmox_esm_Antarctica.nml after ismip7_ant
# (scripts/ismip7/params.sh): what 1pctCO2 changes beyond the ISMIP7 setup.
pct_ant=(
    "transient.tstep_const=2020"            # [yr] time of the "const" method
    "transient.time_init=2020.0"            # [yr] 1pctCO2 start (run year 0 = 1xCO2)
    "transient.time_end=2160.0"             # [yr] 1pctCO2 end (140 yr; ~4xCO2 near yr 140)
    "transient.time_esm_ref==2020,2021"     # [yr] self-reference: dts = tas(t) - mean(tas over 2020-2021)
    "transient.time_hist==2018,2019"        # [yr] no historical period: kept below time_init (use_hist=False)
    "transient.time_proj==2020,2160"        # [yr] projection window (full run)
    "comps.surface_method=chion"            # the ESM gives no smb
    # present-day smb and T_srf: RACMO2.3 (VW23)
    "yelmo_data.pd_tsrf_path=ice_data/{domain}/{grid_name}/RACMO2.3/{grid_name}_RACMO23-VW23.nc"
    "yelmo_data.pd_smb_path=ice_data/{domain}/{grid_name}/RACMO2.3/{grid_name}_RACMO23-VW23.nc"
    "esm.par_file=input/esm/esm_ant_1pctCO2.nml"
    "esm.experiment=1pctCO2"                # the spin-up sets ctrl
    "esm.esm_name=MPI-ESM1-2-LR"            # {gcm} of the forcing files
)
