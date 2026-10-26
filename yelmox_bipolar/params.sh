# yelmox_bipolar coupling setups on yelmox_bipolar/yelmox_bipolar_Bipolar.nml,
# as runme -p arrays (docs: yelmox-bipolar.md, "Coupling setups"). Sourced by
# run_tests.sh of this folder, which runs each of them. Use:
#
#     source yelmox_bipolar/params.sh
#     runme -rs -e bipolar -n yelmox_bipolar/yelmox_bipolar_Bipolar.nml \
#         -o output/bip-full -p "${bip_full[@]}"
#
# The par file itself is the uncoupled setup. The exchanges with the ocean box
# model (OBM, [ctrl] active_obm = True):
#   atmosphere -> OBM  [ctrl] atm2obm: gmt -> the box air temperatures and
#                      vapour fluxes ([obm_forcing])
#   ice sheet  -> OBM  [ctrl] couple_fwf_north/south: the freshwater flux of a
#                      domain into its box (needs the hydrographic masks
#                      [ctrl] hydro_mask_north/south)
#   OBM -> ice sheet   the term obm_north/obm_south of a domain: its ocean
#                      temperature from its box, in place of the term ocean;
#                      the coefficient of the same name in [forcing] coefs
# The OBM does not act on the atmosphere of the domains: gmt alone drives it.
# Other combinations follow from these, e.g. bip_full without
# ctrl.couple_fwf_south=True: Greenland <-> OBM -> Antarctica.

# 1. Uncoupled: two independent ice sheets, no box model (the par file).
bip_uncoupled=()

# 2. The box model alone, forced by gmt.
bip_obm_only=(
    "ctrl.active_north=False"
    "ctrl.active_south=False"
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
)

# 3. OBM -> both ice sheets (one-way): the ocean of each domain from its box.
bip_obm_to_ice=(
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
    "forcing.coefs==gmt,obm_north,obm_south"
    "forcing_north.terms==regional_north,obm_north"
    "forcing_south.terms==regional_south,obm_south"
)

# 4. Greenland -> OBM (one-way): the freshwater flux of Greenland into the
# northern box; the ocean of Greenland as in the par file (term ocean).
bip_gis_to_obm=(
    "ctrl.active_south=False"
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
    "ctrl.couple_fwf_north=True"
)

# 5. Greenland <-> OBM.
bip_gis_obm=(
    "ctrl.active_south=False"
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
    "ctrl.couple_fwf_north=True"
    "forcing.coefs==gmt,obm_north"
    "forcing_north.terms==regional_north,obm_north"
)

# 6. Antarctica <-> OBM.
bip_ais_obm=(
    "ctrl.active_north=False"
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
    "ctrl.couple_fwf_south=True"
    "forcing.coefs==gmt,obm_south"
    "forcing_south.terms==regional_south,obm_south"
)

# 7. Greenland <-> OBM <-> Antarctica (fully coupled).
bip_full=(
    "ctrl.active_obm=True"
    "ctrl.atm2obm=True"
    "ctrl.couple_fwf_north=True"
    "ctrl.couple_fwf_south=True"
    "forcing.coefs==gmt,obm_north,obm_south"
    "forcing_north.terms==regional_north,obm_north"
    "forcing_south.terms==regional_south,obm_south"
)

# Modifier, after a setup with obm_north/obm_south: the ocean as the reference
# plus the box anomaly (theta - theta_init), instead of the box temperature at
# every depth. E.g. -p "${bip_full[@]}" "${bip_ocean_anomaly[@]}".
bip_ocean_anomaly=(
    "coef_obm_north.anomaly=True"
    "term_obm_north.pattern=uniform"
    "coef_obm_south.anomaly=True"
    "term_obm_south.pattern=uniform"
)

# The setups run by run_tests.sh, in order.
bip_setups=(
    bip_uncoupled
    bip_obm_only
    bip_obm_to_ice
    bip_gis_to_obm
    bip_gis_obm
    bip_ais_obm
    bip_full
)
