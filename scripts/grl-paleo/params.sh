# Greenland paleo forcing on the base yelmox/yelmox_GreenlandPaleo.nml, as runme
# -p arrays. Sourced by the run_*.sh scripts of this folder and by
# scripts/docs-figures/run.sh. The named groups (series_*, term_*, snapshot_*)
# come from input/forcing_db.nml (docs/_forcing-terms.md). run_spinup.sh takes
# any of them by name (grlpal_<forcing>).

# Reference climate only (the present-day snapshot of the base par file).
grlpal_pd=()

# A: last glacial period (run_lgp.sh; run_spinup.sh lgp). PMIP3 piControl and
# LGM blended by the glacial index input/alpha_combined_125kyr_interp.dat.
grlpal_lgp=(
    "forcing.series==alpha,dTo"         # alpha: the glacial index ([series_alpha] of the base par file); dTo: 0.25 x the mean air anomaly
    "forcing.terms==lgm,ocean"          # lgm: lgm - picontrol (PMIP3) weighted by alpha, on tas and pr;
                                        # ocean: the profiles plus dTo
    "orbit.method=forcing_time"                # transient insolation (orbit at the forcing time)
)

# B: deglaciation after Tabone et al. (2024) (run_spinup.sh recon, run_deglac.sh).
# Anomalies of the Greenland reconstructions on the present-day reference:
# temperature of Buizert et al. (2018, axis b18: reconstruction - its 1981-2010
# baseline), precipitation of Badgeley et al. (2020, axis b20: the fraction
# itself), tapered outside the GRL box. Swap the fields of term_b18 and
# term_b20 (-p term_b18.fields=pr term_b20.fields=tas) for the other combinations.
grlpal_recon=(
    "forcing.series==dTo"                   # dTo: 0.25 x the mean air anomaly (derived)
    "forcing.terms==b18,b20,ocean"          # b18: Buizert temperature anomaly in time (tas), in full;
                                            # b20: Badgeley precipitation ratio in time (pr), in full;
                                            # ocean: the profiles plus dTo
    "taper.active=True"                     # towards the reference outside the GRL box:
    "taper.weight=taper_b18_b20"            #   w = 1 in the reconstructions' coverage, 0 beyond 100 km
    "taper.outside=mean"                    #   where w < 1: ref + the mean anomaly of the cells with w = 1
    "orbit.method=forcing_time"                    # insolation at the forcing time
)
