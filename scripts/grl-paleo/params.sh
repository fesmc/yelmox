# Greenland paleo forcing on the base yelmox/yelmox_GreenlandPaleo.nml, as runme
# -p arrays. Sourced by the run_*.sh scripts of this folder and by
# scripts/docs-figures/run.sh. The named groups (idx_*, snap_*, snap_axis_*) come
# from input/forcing_db.nml.

# A: last glacial period (run_lgp.sh). PMIP3 piControl and LGM blended by the
# glacial index input/alpha_combined_125kyr_interp.dat.
grlpal_lgp=(
    "indices.names==a"                  # glacial index [idx_a] (series)
    "snap.snapshots==pd,picontrol,lgm"  # present day plus PMIP3 piControl and LGM
    "snap.axes==a"                      # axis picontrol -> lgm, position from index a
    "snap_field_tas.axes=a"             # temperature on the axis a
    "snap_field_pr.blend=ratio"         # precipitation as a ratio
    "snap_field_pr.axes=a"              # ... on the axis a
    "ocn_fraction.active=True"          # ocean anomaly from the air temperature anomaly
    "orbit.method=clock"                # transient insolation (orbit at the forcing clock)
)

# B: deglaciation after Tabone et al. (2024) (run_lgm_spinup.sh, run_deglac.sh).
# Anomalies of the Greenland reconstructions on the present-day reference:
# temperature of Buizert et al. (2018, axis b18: reconstruction - its 1981-2010
# baseline), precipitation of Badgeley et al. (2020, axis b20: the fraction
# itself), tapered outside the GRL box. Swap the axes of snap_field_tas and
# snap_field_pr for the other combinations.
grlpal_recon=(
    "indices.names==one"                    # constant 1: the full anomaly ([idx_one])
    "snap.snapshots==pd,b18_base,b18,b20"   # present day plus the reconstructions
    "snap.axes==b18,b20"                    # axes b18 (b18_base -> b18) and b20 ("" -> b20)
    "snap_field_tas.axes=b18"               # temperature: b18 (Buizert) or b20 (Badgeley)
    "snap_field_pr.blend=ratio"             # precipitation as a ratio
    "snap_field_pr.axes=b20"                # ... b20 (Badgeley) or b18 (Buizert)
    "taper.active=True"                     # towards the reference outside the GRL box:
    "taper.var_defs=input/clim_db.nml"      #   database of the weight
    "taper.weight=taper_b18_b20"            #   w = 1 in the reconstructions' coverage, 0 beyond 100 km
    "taper.outside=mean"                    #   where w < 1: ref + the mean anomaly of the cells with w = 1
    "ocn_fraction.active=True"              # ocean anomaly from the air temperature anomaly
    "orbit.method=clock"                    # insolation at the forcing clock
)
