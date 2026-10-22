module obm_coupling
    ! Bipolar ocean coupling: exchanges scalars between the kryos_domains (Yelmo
    ! + climate) and the shared Ocean Box Model (OBM). This is a bridge module --
    ! it sits above both the Kryos modules (kryos_domain) and the libs/obm ocean
    ! box model. It is only pertinent to yelmox_bipolar, so it lives here next
    ! to its driver.
    !
    ! The whole obm stack is single precision (obm_defs preci = kind(1.0)), which
    ! matches Yelmo's wp (= sp), so obm and Yelmo/climate fields are exchanged
    ! without kind conversion.
    !
    ! The three exchanges:
    !   atm2obm  : global mean air temperature anomaly ([obm_forcing], a
    !              coefficient of the driver's registry) -> obm box air temperatures/vapor
    !   ism2obm  : ice-sheet freshwater flux (calc_fwf)  -> obm%fn / obm%fs
    !   obm2ism  : obm box ocean temperature -> the domain's ocean temperature
    !              to_ann (marine-shelf forcing), through the coefficients
    !              derived from the boxes (obm_north, obm_south) and the domains' ocean
    !              terms on them (climate_terms), applied here
    !
    ! The module owns its whole configuration (obm_coupling_ctl: exchange
    ! switches, the atmospheric forcing, freshwater-flux masks, nautilus
    ! hysteresis forcing), loaded by obm_ctl_load from [ctrl] and
    ! [obm_forcing]. obm2ism is on when a coefficient is derived from a box.
    ! The driver calls obm_ctl_load, then coupling_obm2ism once the box model
    ! is initialized (cold start or restart) and before the domains, so their
    ! startup climate sees the box, then obm_coupling_init once the domains are
    ! initialized, then obm_exchange once per step after the domains' climate
    ! steps.

    use nml,           only : nml_read, nml_validate, nml_replace
    use ncio,          only : nc_read
    use precision,     only : dp
    use yelmo,         only : wp
    use kryos,         only : kryos_domain
    use kryos_forcing, only : coefs_class, coefs_has, coefs_value, coefs_uses_derived, &
                              coefs_set_obm
    use yelmox_climate, only : climate_apply_ocean
    use yelmox_nml,    only : def_file => yelmox_defaults_file
    use obm_defs,      only : obm_class
    use ice2ocean,     only : calc_fwf
    use ocean2ice,     only : obm_temperature_kelvin

    implicit none
    private

    type obm_coupling_ctl
        ! Which exchanges are active ([ctrl]; obm2ism: a coefficient derived from a box).
        logical :: active_obm = .false.
        logical :: ism2obm    = .false.
        logical :: atm2obm    = .false.
        character(len=512) :: obm_name = "none"   ! OBM parameter block ("nautilus"/"stommel")
        ! atm2obm: the global mean air temperature anomaly dT_glob = gain * coef
        ! ([obm_forcing], a coefficient).
        character(len=64)  :: forcing_coef  = ""
        real(wp)           :: forcing_gain  = 1.0_wp
        ! ism2obm: freshwater-flux definition + per-hemisphere hydrographic masks.
        logical            :: couple_fwf_north = .false.
        logical            :: couple_fwf_south = .false.
        character(len=512) :: fwf_definition   = "dVdt"
        character(len=512) :: hydro_mask_north_path = ""
        character(len=512) :: hydro_mask_south_path = ""
        real(wp), allocatable :: hydro_mask_north(:,:), hydro_mask_south(:,:)
        ! Hysteresis forcing (nautilus only, from the [nautilus] group).
        logical            :: hyster_on = .false.
        character(len=512) :: hyster_forcing = "", hyster_forcing_method = ""
        real(wp)           :: hyster_rate = 0.0_wp
        real(dp)           :: hyster_positive_branch_time = 0.0_dp
    end type obm_coupling_ctl

    public :: obm_coupling_ctl, obm_ctl_load, obm_coupling_init, obm_exchange
    public :: coupling_atm2obm, coupling_ism2obm, coupling_obm2ism
    public :: update_bipolar_hyster_forcing

contains

    subroutine obm_ctl_load(oc, path_par)
        ! Load the ocean-coupling switches ([ctrl]) and, when needed, the
        ! atmospheric forcing ([obm_forcing]; unset keys from
        ! input/yelmox_defaults.nml), the freshwater-flux mask paths and the
        ! nautilus hysteresis-forcing parameters. The masks themselves are read
        ! by obm_coupling_init (they need the initialized domains for their
        ! grid shape).
        type(obm_coupling_ctl), intent(inout) :: oc
        character(len=*),       intent(in)    :: path_par

        call nml_read(path_par, "ctrl", "active_obm", oc%active_obm)
        call nml_read(path_par, "ctrl", "ism2obm",    oc%ism2obm)
        call nml_read(path_par, "ctrl", "atm2obm",    oc%atm2obm)
        call nml_read(path_par, "ctrl", "obm_name",   oc%obm_name)

        if (oc%atm2obm) then
            call nml_validate(path_par, def_file, "obm_forcing")
            call nml_read(path_par, "obm_forcing", "coef",  oc%forcing_coef, defaults_file=def_file)
            call nml_read(path_par, "obm_forcing", "gain",  oc%forcing_gain,  defaults_file=def_file)
        end if

        if (oc%ism2obm) then
            call nml_read(path_par, "ctrl", "couple_fwf_north", oc%couple_fwf_north)
            call nml_read(path_par, "ctrl", "couple_fwf_south", oc%couple_fwf_south)
            call nml_read(path_par, "ctrl", "fwf_definition",   oc%fwf_definition)
            if (oc%couple_fwf_north) &
                call nml_read(path_par, "ctrl", "hydro_mask_north", oc%hydro_mask_north_path)
            if (oc%couple_fwf_south) &
                call nml_read(path_par, "ctrl", "hydro_mask_south", oc%hydro_mask_south_path)
        end if

        if (trim(oc%obm_name) == "nautilus") then
            call nml_read(path_par, oc%obm_name, "hyster_on",                   oc%hyster_on)
            call nml_read(path_par, oc%obm_name, "hyster_forcing",              oc%hyster_forcing)
            call nml_read(path_par, oc%obm_name, "hyster_forcing_method",       oc%hyster_forcing_method)
            call nml_read(path_par, oc%obm_name, "hyster_rate",                 oc%hyster_rate)
            call nml_read(path_par, oc%obm_name, "hyster_positive_branch_time", oc%hyster_positive_branch_time)
        end if
    end subroutine obm_ctl_load

    subroutine obm_coupling_init(oc, coefs, dom_north, dom_south, active_north, active_south)
        ! Check the exchanges against what the run holds, then load the
        ! hydrographic masks (Yelmo grid) restricting the freshwater flux per
        ! hemisphere. Every exchange needs the box model (active_obm); atm2obm
        ! the coefficient of [obm_forcing] in the registry; ism2obm an ice-sheet
        ! domain. Only touches a domain when it is active, so an inactive domain
        ! may be uninitialized. Call after the coefficients and the domains.
        type(obm_coupling_ctl), intent(inout) :: oc
        type(coefs_class),      intent(in)    :: coefs
        type(kryos_domain),     intent(in)    :: dom_north, dom_south
        logical,                intent(in)    :: active_north, active_south

        logical :: obm2ism

        ! obm2ism: a coefficient derived from a box (obm_north, obm_south), used by
        ! the ocean terms of a domain.
        obm2ism = coefs_uses_derived(coefs, "obm_north") .or. coefs_uses_derived(coefs, "obm_south")

        if (.not. oc%active_obm .and. (oc%atm2obm .or. oc%ism2obm .or. obm2ism)) then
            write(*,*) "obm_coupling_init:: error: atm2obm, ism2obm and the coefficients derived &
                       &from the boxes (obm_north, obm_south) need the ocean box model &
                       &([ctrl] active_obm = True)."
            error stop 1
        end if
        if (oc%atm2obm .and. .not. coefs_has(coefs, oc%forcing_coef)) then
            write(*,*) "obm_coupling_init:: error: [obm_forcing] coef = '"// &
                       trim(oc%forcing_coef)//"' is not in [forcing] coefs (needed by atm2obm)."
            error stop 1
        end if
        if (oc%ism2obm .and. .not. (active_north .or. active_south)) then
            write(*,*) "obm_coupling_init:: error: ism2obm needs active_north or active_south."
            error stop 1
        end if

        if (oc%atm2obm) &
            write(*,*) "obm_coupling_init:: atm2obm from coefficient "//trim(oc%forcing_coef)// &
                       ", gain =", oc%forcing_gain

        if (.not. oc%ism2obm) return

        ! The masks are on the Yelmo grid of their domain: {domain}/{grid_name}
        ! of the paths are its name and grid_ice.
        if (oc%couple_fwf_north .and. active_north) then
            call nml_replace(oc%hydro_mask_north_path, "{domain}",    trim(dom_north%ctl%domain))
            call nml_replace(oc%hydro_mask_north_path, "{grid_name}", trim(dom_north%ctl%grid_ice))
            allocate(oc%hydro_mask_north(dom_north%yelmo%grd%G%nx, dom_north%yelmo%grd%G%ny))
            call nc_read(oc%hydro_mask_north_path, "mask", oc%hydro_mask_north)
        end if
        if (oc%couple_fwf_south .and. active_south) then
            call nml_replace(oc%hydro_mask_south_path, "{domain}",    trim(dom_south%ctl%domain))
            call nml_replace(oc%hydro_mask_south_path, "{grid_name}", trim(dom_south%ctl%grid_ice))
            allocate(oc%hydro_mask_south(dom_south%yelmo%grd%G%nx, dom_south%yelmo%grd%G%ny))
            call nc_read(oc%hydro_mask_south_path, "mask", oc%hydro_mask_south)
        end if
    end subroutine obm_coupling_init

    subroutine obm_exchange(oc, obm, coefs, dom_north, dom_south, active_north, active_south, &
                            time, time_init, dtt)
        ! One coupling pass, in the yelmox_bipolar order: atmosphere -> obm,
        ! ice-sheet freshwater flux -> obm, hysteresis forcing (nautilus), then
        ! obm ocean temperature -> ice sheets (read by the marine-shelf step).
        ! Call after the domains' dynamics/climate have advanced and before
        ! their marine-shelf steps. coefs: the driver's registry, advanced to
        ! this step.
        type(obm_coupling_ctl), intent(inout) :: oc
        type(obm_class),        intent(inout) :: obm
        type(coefs_class),      intent(inout) :: coefs
        type(kryos_domain),     intent(inout) :: dom_north, dom_south
        logical,                intent(in)    :: active_north, active_south
        real(dp),               intent(in)    :: time, time_init, dtt

        if (oc%atm2obm) &
            call coupling_atm2obm(obm, oc%forcing_gain*coefs_value(coefs, oc%forcing_coef))

        if (oc%ism2obm) then
            if (oc%couple_fwf_north .and. active_north) &
                call coupling_ism2obm(dom_north, obm, oc%hydro_mask_north, "north", oc%fwf_definition)
            if (oc%couple_fwf_south .and. active_south) &
                call coupling_ism2obm(dom_south, obm, oc%hydro_mask_south, "south", oc%fwf_definition)
        end if

        if (trim(oc%obm_name) == "nautilus" .and. oc%hyster_on) &
            call update_bipolar_hyster_forcing(time, time_init, obm, dtt, &
                    oc%hyster_positive_branch_time, oc%hyster_rate, &
                    oc%hyster_forcing, oc%hyster_forcing_method)

        ! obm2ism: the coefficients derived from the boxes, then the ocean terms of
        ! each domain re-applied (idempotent; a no-op without ocean terms).
        call coupling_obm2ism(obm, oc%obm_name, coefs)
        if (active_north) call climate_apply_ocean(dom_north%cl, dom_north%clim, coefs)
        if (active_south) call climate_apply_ocean(dom_south%cl, dom_south%clim, coefs)
    end subroutine obm_exchange

    subroutine coupling_atm2obm(obm, dT_glob)
        ! Atmosphere -> OBM: the box-model air temperatures of every box from the
        ! global mean air temperature anomaly dT_glob [K] with the regional
        ! amplification factors (Zickfeld et al., 2004), theta_i = theta_i_init
        ! + ampl_i * dT_glob (north, tropics, south), and the vapor fluxes,
        ! phin (tropics -> north) and phit (south -> tropics), linear in dT_glob
        ! by the hydrological sensitivities.
        type(obm_class), intent(inout) :: obm
        real(wp),        intent(in)    :: dT_glob

        obm%thetan = obm%par%thetan_init + obm%par%thermal_ampl_north  *dT_glob
        obm%thetat = obm%par%thetat_init + obm%par%thermal_ampl_tropics*dT_glob
        obm%thetas = obm%par%thetas_init + obm%par%thermal_ampl_south  *dT_glob
        obm%phin   = obm%par%phin_init   + obm%par%hn*obm%par%pnh*dT_glob
        obm%phit   = obm%par%phit_init   + obm%par%hs*obm%par%psh*dT_glob
    end subroutine coupling_atm2obm

    subroutine coupling_ism2obm(dom, obm, mask, hemisphere, fwf_def)
        ! Ice sheet -> OBM: freshwater flux from the domain's ice mass balance
        ! (calc_fwf, on the Yelmo grid, restricted by the hydrographic mask) into
        ! the box model's northern (fn) or southern (fs) input flux.
        ! dom is intent(inout) because calc_fwf takes its mass-balance fields as
        ! non-intent (modifiable) allocatable dummies.
        type(kryos_domain),      intent(inout) :: dom
        type(obm_class),       intent(inout) :: obm
        real(wp), allocatable, intent(inout) :: mask(:,:)
        character(len=*),      intent(in)    :: hemisphere
        character(len=*),      intent(in)    :: fwf_def

        real(wp)           :: fwf
        character(len=512) :: fdef

        ! calc_fwf reassigns its fwf_def dummy internally, so pass a local copy.
        fdef = fwf_def
        fwf  = calc_fwf(dom%yelmo%bnd%c%rho_w, dom%yelmo%bnd%c%rho_ice, &
                        dom%yelmo%bnd%c%sec_year, &
                        dom%yelmo%tpo%now%mb_net, dom%yelmo%tpo%now%smb, &
                        dom%yelmo%tpo%now%bmb,    dom%yelmo%tpo%now%cmb, &
                        dom%yelmo%tpo%now%H_ice,  dom%yelmo%tpo%now%dHidt, &
                        dom%yelmo%tpo%now%f_grnd, &
                        dom%yelmo%tpo%par%dx, dom%yelmo%tpo%par%dy, &
                        mask, hemisphere, fdef)

        select case(trim(hemisphere))
            case("north"); obm%fn = fwf
            case("south"); obm%fs = fwf
        end select
    end subroutine coupling_ism2obm

    subroutine coupling_obm2ism(obm, obm_name, coefs)
        ! OBM -> ice sheets: the box-model ocean temperatures of the
        ! hemispheres (northern box tn, southern box ts) and their initial
        ! values, in K, into the coefficients derived from them (obm_north,
        ! obm_south: the temperature, or its anomaly), which the ocean terms of
        ! the domains read (climate_apply_ocean) to set the ocean temperature
        ! that marine_shelf then reads as its forcing.
        type(obm_class),     intent(in)    :: obm
        character(len=*),    intent(in)    :: obm_name
        type(coefs_class),   intent(inout) :: coefs

        call coefs_set_obm(coefs, "obm_north", obm_temperature_kelvin(obm%tn, obm_name), &
                                               obm_temperature_kelvin(obm%par%tn_init, obm_name))
        call coefs_set_obm(coefs, "obm_south", obm_temperature_kelvin(obm%ts, obm_name), &
                                               obm_temperature_kelvin(obm%par%ts_init, obm_name))
    end subroutine coupling_obm2ism

    subroutine update_bipolar_hyster_forcing(t, t0, obm, dt, branch_time_thr, &
                                             rate, forcing, forc_method)
        ! Hysteresis forcing for the nautilus box model: nudge one obm control
        ! (phit/phin/fs/fn, or fn+fs) along a prescribed path.
        real(dp),          intent(in)    :: t, t0
        type(obm_class),   intent(inout) :: obm
        real(dp),          intent(in)    :: dt
        real(dp),          intent(in)    :: branch_time_thr
        real(wp),          intent(in)    :: rate
        character(len=*),  intent(in)    :: forcing, forc_method

        real(wp) :: factor

        factor = 0.0_wp   ! unknown forc_method -> no forcing (safer than the original's UB)

        select case(trim(forc_method))
            case("triangular")
                if (t <= branch_time_thr) then
                    factor = rate * dt
                else
                    factor = -1 * rate * dt
                end if
            case("sin")
                if (trim(forcing) .eq. "phit") then
                    factor = -obm%phit + rate * sin(2*3.14159265358979*(t-t0)/branch_time_thr)
                else if (trim(forcing) .eq. "fn") then
                    factor = -obm%fn + rate * sin(2*3.14159265358979*(t-t0)/branch_time_thr)
                end if
            case("noise")
                if (trim(forcing) .eq. "phit") then
                    factor = -obm%phit + r8_normal_ab(0.0_wp, rate)
                else if (trim(forcing) .eq. "fn") then
                    factor = -obm%fn + r8_normal_ab(0.0_wp, rate)
                end if
            case("linear")
                factor = rate * dt
        end select

        select case(trim(forcing))
            case("phit");        obm%phit = obm%phit + factor
            case("phin");        obm%phin = obm%phin + factor
            case("fs");          obm%fs   = obm%fs   + factor
            case("fn");          obm%fn   = obm%fn   + factor
            case("global_melt")
                obm%fn = obm%fn + factor
                obm%fs = obm%fs + factor
        end select
    end subroutine update_bipolar_hyster_forcing

    ! ----- private helpers -----


    function r8_normal_ab(a, b) result(val)
        ! Sample of a normal PDF with mean a, standard deviation b (Box-Muller).
        ! After John Burkardt (MIT license).
        real(wp), intent(in) :: a, b
        real(wp) :: val

        integer, parameter :: rk = kind(1.0D+00)
        real(kind=rk) :: r1, r2, x
        real(kind=rk), parameter :: r8_pi = 3.141592653589793D+00

        call random_number(harvest=r1)
        call random_number(harvest=r2)
        x = sqrt(-2.0D+00 * log(r1)) * cos(2.0D+00 * r8_pi * r2)
        val = a + b * real(x, wp)
    end function r8_normal_ab

end module obm_coupling
