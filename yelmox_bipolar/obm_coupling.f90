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
    !   atm2obm  : air temperature anomalies ([obm_forcing], coefficients of the
    !              driver's registry: one per box, and the global one) -> obm
    !              box air temperatures/vapor
    !   ism2obm  : ice-sheet freshwater flux (calc_fwf)  -> obm%fn / obm%fs,
    !              per hemisphere ([ctrl] couple_fwf_north, couple_fwf_south)
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
                              coefs_set_obm, coefs_derived, derived_spec, coefs_varies
    use yelmox_climate, only : climate_apply_ocean
    use yelmox_nml,    only : def_file => yelmox_defaults_file
    use obm_defs,      only : obm_class
    use ice2ocean,     only : calc_fwf
    use ocean2ice,     only : obm_temperature_kelvin

    implicit none
    private

    type obm_coupling_ctl
        ! Which exchanges are active ([ctrl]; ism2obm: couple_fwf_north/south below;
        ! obm2ism: a coefficient derived from a box).
        logical :: active_obm = .false.
        logical :: atm2obm    = .false.
        character(len=512) :: obm_name = "none"   ! OBM parameter block ("nautilus")
        ! atm2obm: the coefficients of [obm_forcing], the air temperature
        ! anomaly of each box and the global one (the vapour fluxes) [K].
        character(len=64)  :: coef_north = "", coef_tropics = "", coef_south = "", coef_global = ""
        ! ism2obm: the freshwater flux of each hemisphere into its box, its
        ! definition and the hydrographic masks.
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

    ! The keys of [obm_forcing], in the order of obm_forcing_coefs.
    character(len=12), parameter :: OBM_FORCING_KEYS(4) = &
        [character(len=12) :: "coef_north", "coef_tropics", "coef_south", "coef"]

    public :: obm_coupling_ctl, obm_ctl_load, obm_coupling_init, obm_exchange, obm_forcing_coefs
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

        call nml_read(path_par, "ctrl", "active_obm",       oc%active_obm)
        call nml_read(path_par, "ctrl", "atm2obm",          oc%atm2obm)
        call nml_read(path_par, "ctrl", "couple_fwf_north", oc%couple_fwf_north)
        call nml_read(path_par, "ctrl", "couple_fwf_south", oc%couple_fwf_south)
        call nml_read(path_par, "ctrl", "obm_name",         oc%obm_name)

        if (oc%atm2obm) then
            call nml_validate(path_par, def_file, "obm_forcing")
            call nml_read(path_par, "obm_forcing", "coef_north",   oc%coef_north,   defaults_file=def_file)
            call nml_read(path_par, "obm_forcing", "coef_tropics", oc%coef_tropics, defaults_file=def_file)
            call nml_read(path_par, "obm_forcing", "coef_south",   oc%coef_south,   defaults_file=def_file)
            call nml_read(path_par, "obm_forcing", "coef",         oc%coef_global,  defaults_file=def_file)
        end if

        if (oc%couple_fwf_north .or. oc%couple_fwf_south) then
            call nml_read(path_par, "ctrl", "fwf_definition",   oc%fwf_definition)
            select case(trim(oc%fwf_definition))
                case("dVdt", "mb", "ablation")
                case default
                    write(*,*) "obm_ctl_load:: error: [ctrl] fwf_definition must be dVdt, mb &
                               &or ablation; got '"//trim(oc%fwf_definition)//"'."
                    error stop 1
            end select
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
        ! the coefficients of [obm_forcing] in the registry; the freshwater flux
        ! of a hemisphere its domain. The ocean terms of the domains against the
        ! shared coefficients (check_domain_ocean, and an air-anomaly coefficient
        ! set by one domain only). Only touches a domain when it is active, so
        ! an inactive domain may be uninitialized. Call after the coefficients
        ! and the domains.
        type(obm_coupling_ctl), intent(inout) :: oc
        type(coefs_class),      intent(in)    :: coefs
        type(kryos_domain),     intent(in)    :: dom_north, dom_south
        logical,                intent(in)    :: active_north, active_south

        logical :: obm2ism
        character(len=64), allocatable :: names(:)
        integer :: k

        ! obm2ism: a coefficient derived from a box (obm_north, obm_south), used by
        ! the ocean terms of a domain.
        obm2ism = coefs_uses_derived(coefs, "obm_north") .or. coefs_uses_derived(coefs, "obm_south")

        if (.not. oc%active_obm .and. (oc%atm2obm .or. oc%couple_fwf_north .or. &
                                       oc%couple_fwf_south .or. obm2ism)) then
            write(*,*) "obm_coupling_init:: error: atm2obm, couple_fwf_north/south and the &
                       &coefficients derived from the boxes (obm_north, obm_south) need the &
                       &ocean box model ([ctrl] active_obm = True)."
            error stop 1
        end if
        names = obm_forcing_coefs(oc)
        do k = 1, size(names)
            if (coefs_has(coefs, names(k))) cycle
            write(*,*) "obm_coupling_init:: error: [obm_forcing] "//trim(OBM_FORCING_KEYS(k))// &
                       " = '"//trim(names(k))//"' is not in [forcing] coefs (needed by atm2obm)."
            error stop 1
        end do
        if (oc%couple_fwf_north .and. .not. active_north) then
            write(*,*) "obm_coupling_init:: error: couple_fwf_north needs active_north."
            error stop 1
        end if
        if (oc%couple_fwf_south .and. .not. active_south) then
            write(*,*) "obm_coupling_init:: error: couple_fwf_south needs active_south."
            error stop 1
        end if

        if (active_north) call check_domain_ocean(dom_north, coefs, "north")
        if (active_south) call check_domain_ocean(dom_south, coefs, "south")

        ! The box model sees the forcing of the ice sheets: with an atmosphere
        ! of a domain that changes in time, atm2obm forces the boxes too.
        if (oc%active_obm .and. .not. oc%atm2obm) then
            if (active_north) call check_box_forced(dom_north, coefs, "north")
            if (active_south) call check_box_forced(dom_south, coefs, "south")
        end if

        ! A coefficient derived from the air anomaly is set by the climate of
        ! the domain whose terms use it. The coefficients are shared and
        ! obm_exchange re-applies the ocean terms of both domains after both
        ! climate updates, so a coefficient set by both would give the first
        ! domain the anomaly of the second.
        if (active_north .and. active_south .and. allocated(dom_north%cl%dair) .and. &
            allocated(dom_south%cl%dair)) then
            do k = 1, size(dom_north%cl%dair)
                if (any(dom_south%cl%dair%name == dom_north%cl%dair(k)%name)) then
                    write(*,*) "obm_coupling_init:: error: the coefficient "// &
                               trim(dom_north%cl%dair(k)%name)//" (derived from the air &
                               &anomaly) is used by the terms of both domains; the &
                               &coefficients are shared, so each domain needs a coefficient &
                               &of its own (e.g. [coef_dTo_south] as [coef_dTo])."
                    error stop 1
                end if
            end do
        end if

        if (oc%atm2obm) &
            write(*,*) "obm_coupling_init:: atm2obm from the coefficients "//trim(oc%coef_north)// &
                       ", "//trim(oc%coef_tropics)//", "//trim(oc%coef_south)//" (boxes), "// &
                       trim(oc%coef_global)//" (vapour fluxes)"

        ! The masks are on the Yelmo grid of their domain: {domain}/{grid_name}
        ! of the paths are its name and grid_ice.
        if (oc%couple_fwf_north) then
            call nml_replace(oc%hydro_mask_north_path, "{domain}",    trim(dom_north%ctl%domain))
            call nml_replace(oc%hydro_mask_north_path, "{grid_name}", trim(dom_north%ctl%grid_ice))
            allocate(oc%hydro_mask_north(dom_north%yelmo%grd%G%nx, dom_north%yelmo%grd%G%ny))
            call nc_read(oc%hydro_mask_north_path, "mask", oc%hydro_mask_north)
        end if
        if (oc%couple_fwf_south) then
            call nml_replace(oc%hydro_mask_south_path, "{domain}",    trim(dom_south%ctl%domain))
            call nml_replace(oc%hydro_mask_south_path, "{grid_name}", trim(dom_south%ctl%grid_ice))
            allocate(oc%hydro_mask_south(dom_south%yelmo%grd%G%nx, dom_south%yelmo%grd%G%ny))
            call nc_read(oc%hydro_mask_south_path, "mask", oc%hydro_mask_south)
        end if
    end subroutine obm_coupling_init

    subroutine check_domain_ocean(dom, coefs, hemisphere)
        ! The ocean terms of a domain against the box model: a coefficient
        ! derived from a box is that of the domain's hemisphere (obm_north for
        ! the north), and an ocean from a box takes no anomaly from the air on
        ! top (the box follows the air already, through atm2obm).
        type(kryos_domain), intent(in) :: dom
        type(coefs_class),  intent(in) :: coefs
        character(len=*),   intent(in) :: hemisphere    ! north | south

        type(derived_spec) :: drv
        character(len=64)  :: box_term, air_term
        integer :: k

        box_term = ""
        air_term = ""
        do k = 1, dom%cl%terms%n
            associate(t => dom%cl%terms%term(k))
            if (.not. (t%on_to .and. t%derived)) cycle
            drv = coefs_derived(coefs, t%coef)
            select case(trim(drv%from))
                case("obm_north", "obm_south")
                    if (trim(drv%from) /= "obm_"//trim(hemisphere)) then
                        write(*,*) "obm_coupling_init:: error: ["//trim(t%group)//"] coef "// &
                                   trim(t%coef)//" is derived from the box "//trim(drv%from)// &
                                   "; the "//trim(hemisphere)//" domain takes its ocean from &
                                   &obm_"//trim(hemisphere)//"."
                        error stop 1
                    end if
                    box_term = t%group
                case("air_anomaly")
                    air_term = t%group
            end select
            end associate
        end do

        if (len_trim(box_term) > 0 .and. len_trim(air_term) > 0) then
            write(*,*) "obm_coupling_init:: error: the "//trim(hemisphere)//" domain has both ["// &
                       trim(box_term)//"] (its ocean from the box) and ["//trim(air_term)// &
                       "] (an ocean anomaly from the air); the box follows the air already, &
                       &use one of them (e.g. obm_"//trim(hemisphere)//" in place of ocean)."
            error stop 1
        end if
    end subroutine check_domain_ocean

    subroutine check_box_forced(dom, coefs, hemisphere)
        ! Without atm2obm the box model keeps its initial air temperatures: an
        ! error if a term of the domain changes its atmosphere in time.
        type(kryos_domain), intent(in) :: dom
        type(coefs_class),  intent(in) :: coefs
        character(len=*),   intent(in) :: hemisphere    ! north | south
        integer :: k
        do k = 1, dom%cl%terms%n
            associate(t => dom%cl%terms%term(k))
            if (.not. (t%on_tas .or. t%on_pr) .or. len_trim(t%coef) == 0) cycle
            if (.not. coefs_varies(coefs, t%coef)) cycle
            write(*,*) "obm_coupling_init:: error: ["//trim(t%group)//"] of the "//trim(hemisphere)// &
                       " domain changes its atmosphere in time (coef "//trim(t%coef)//"), but the &
                       &box model is not forced by the atmosphere: set [ctrl] atm2obm = True."
            error stop 1
            end associate
        end do
    end subroutine check_box_forced

    function obm_forcing_coefs(oc) result(names)
        ! The coefficients atm2obm reads ([obm_forcing]): the air temperature
        ! anomaly of the north, tropics and south box and the global one; none
        ! without atm2obm.
        type(obm_coupling_ctl), intent(in) :: oc
        character(len=64), allocatable :: names(:)
        if (oc%atm2obm) then
            names = [character(len=64) :: oc%coef_north, oc%coef_tropics, oc%coef_south, oc%coef_global]
        else
            allocate(names(0))
        end if
    end function obm_forcing_coefs

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
            call coupling_atm2obm(obm, coefs_value(coefs, oc%coef_north), coefs_value(coefs, oc%coef_tropics), &
                                  coefs_value(coefs, oc%coef_south), coefs_value(coefs, oc%coef_global))

        if (oc%couple_fwf_north) &
            call coupling_ism2obm(dom_north, obm, oc%hydro_mask_north, "north", oc%fwf_definition)
        if (oc%couple_fwf_south) &
            call coupling_ism2obm(dom_south, obm, oc%hydro_mask_south, "south", oc%fwf_definition)

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

    subroutine coupling_atm2obm(obm, dT_north, dT_tropics, dT_south, dT_glob)
        ! Atmosphere -> OBM: the box-model air temperatures of every box from
        ! its air temperature anomaly [K] (e.g. the global one amplified,
        ! Zickfeld et al., 2004: the scaled coefficients dTa_north, dTa_tropics,
        ! dTa_south), tstar_i = tstar_i_init + dT_i, and the vapor fluxes, phin
        ! (tropics -> north) and phit (south -> tropics), linear in the global
        ! mean air temperature anomaly dT_glob by the hydrological sensitivities.
        type(obm_class), intent(inout) :: obm
        real(wp),        intent(in)    :: dT_north, dT_tropics, dT_south, dT_glob

        obm%tstarn = obm%par%tstarn_init + dT_north
        obm%tstart = obm%par%tstart_init + dT_tropics
        obm%tstars = obm%par%tstars_init + dT_south
        obm%phin   = obm%par%phin_init   + obm%par%hn*obm%par%pnh*dT_glob
        obm%phit   = obm%par%phit_init   + obm%par%hs*obm%par%psh*dT_glob
    end subroutine coupling_atm2obm

    subroutine coupling_ism2obm(dom, obm, mask, hemisphere, fwf_def)
        ! Ice sheet -> OBM: freshwater flux from the domain's ice mass balance
        ! (calc_fwf, on the Yelmo grid, restricted by the hydrographic mask but
        ! for dVdt) into
        ! the box model's northern (fn) or southern (fs) input flux.
        type(kryos_domain),    intent(in)    :: dom
        type(obm_class),       intent(inout) :: obm
        real(wp),              intent(in)    :: mask(:,:)
        character(len=*),      intent(in)    :: hemisphere
        character(len=*),      intent(in)    :: fwf_def

        real(wp) :: fwf

        fwf  = calc_fwf(dom%yelmo%bnd%c%rho_w, dom%yelmo%bnd%c%rho_ice, &
                        dom%yelmo%bnd%c%sec_year, &
                        dom%yelmo%tpo%now%mb_net, dom%yelmo%tpo%now%smb, &
                        dom%yelmo%tpo%now%bmb,    dom%yelmo%tpo%now%cmb, &
                        dom%yelmo%tpo%now%H_ice,  dom%yelmo%tpo%now%dHidt, &
                        dom%yelmo%tpo%now%f_grnd, &
                        dom%yelmo%tpo%par%dx, dom%yelmo%tpo%par%dy, &
                        mask, hemisphere, fwf_def)

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
