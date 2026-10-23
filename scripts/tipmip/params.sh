# TIPMIP setup on the base yelmox/yelmox_esm_Greenland.nml, as a runme -p array.
# Sourced by run_tipmip_greenland.sh and by scripts/docs-figures/run.sh. The
# forcing is the TIPMIP piControl anomalies on top of the MAR + ORAS5 reference
# climatology (input/esm/esm_grl_tipmip.nml); the ESM gives no smb, so the
# surface is chion.

# The ramp (esm-up2p0) on the piControl-branch clock: model years 0..232. The
# stabilisations set their own par_file, time_init and time_end in their step.
tipmip_grl=(
    "transient.tstep_const=0.0"             # [yr] time of the "const" method
    "transient.time_init=0.0"               # [yr] ramp: piControl-branch year 0
    "transient.time_end=232.0"              # [yr] ramp: 232 yr
    "transient.time_ref==1961,1990"         # [yr] reference period of the MAR base climatology
    "transient.time_esm_ref==0,2"           # [yr] ESM reference = ramp start (anomaly ~ 0, ratio ~ 1)
    "transient.time_hist==-1,-1"            # [yr] no historical period: entirely below time_init (use_hist=False)
    "transient.time_proj==0,231"            # [yr] projection = full ramp
    "comps.surface_method=chion"            # the ESM gives no smb
    "esm.par_file=input/esm/esm_grl_tipmip.nml"
    "esm.esm_name=IPSL-CM6-ESMCO2"          # TIPMIP forcing model (also: EC-Earth3-ESM-1, GRL-4KM only)
)
