module obm_coupling
    ! Bipolar ocean coupling: exchanges scalars between an ice_domain (Yelmo +
    ! snapclim) and the shared Ocean Box Model (OBM). Ported from
    ! yelmox_bipolar.f90. This is a bridge module -- it sits above both
    ! yelmox_domain (ice_domain) and the libs/obm ocean box model. It is only
    ! pertinent to the bipolar flavor, so it lives here next to its driver.
    !
    ! The whole obm stack is single precision (obm_defs preci = kind(1.0)), which
    ! matches Yelmo's wp (= sp), so obm and Yelmo/snapclim fields are exchanged
    ! without kind conversion.
    !
    ! The three exchanges (per hemisphere, called from the bipolar driver):
    !   atm2obm  : snapclim air-temperature anomaly -> obm box temperatures/vapor
    !   ism2obm  : ice-sheet freshwater flux (calc_fwf)  -> obm%fn / obm%fs
    !   obm2ism  : obm box ocean temperature -> snapclim%now%to_ann (mshlf forcing)

    use yelmo,         only : wp
    use yelmox_domain, only : ice_domain
    use obm_defs,      only : obm_class
    use ice2ocean,     only : calc_fwf
    use ocean2ice,     only : calc_ocean_temperature_field

    implicit none
    private

    public :: coupling_atm2obm, coupling_ism2obm, coupling_obm2ism
    public :: update_bipolar_hyster_forcing

contains

    subroutine coupling_atm2obm(dom, obm, hemisphere, time)
        ! Atmosphere -> OBM: drive the box-model atmospheric temperatures + vapor
        ! fluxes from this domain's snapclim air-temperature anomaly series.
        ! Hemisphere-specific: north sets thetan/phin, south sets thetas/phit.
        ! Both hemispheres also set the shared tropical box temperature thetat; if
        ! both are active the south value overwrites the north one, exactly as in
        ! yelmox_bipolar (the original flagged this as redundant).
        type(ice_domain), intent(in)    :: dom
        type(obm_class),  intent(inout) :: obm
        character(len=*), intent(in)    :: hemisphere
        real(wp),         intent(in)    :: time

        real(wp) :: at, dTa

        ! Air-temperature anomaly (snapclim series), scaled to a temperature change.
        at  = series_interp(dom%snp%at%time, dom%snp%at%var, time)
        dTa = at * dom%snp%par%dTa_const

        select case(trim(hemisphere))
            case("north")
                obm%thetan = obm%par%thetan_init + dTa
                obm%thetat = obm%par%thetat_init + &
                             dTa*obm%par%thermal_ampl_tropics/obm%par%thermal_ampl_north
                obm%phin   = obm%par%phin_init + &
                             obm%par%hn*obm%par%pnh*dTa/obm%par%thermal_ampl_north
            case("south")
                obm%thetas = obm%par%thetas_init + dTa
                obm%thetat = obm%par%thetat_init + &
                             dTa*obm%par%thermal_ampl_tropics/obm%par%thermal_ampl_south
                obm%phit   = obm%par%phit_init + &
                             obm%par%hs*obm%par%psh*dTa/obm%par%thermal_ampl_south
        end select
    end subroutine coupling_atm2obm

    subroutine coupling_ism2obm(dom, obm, mask, hemisphere, fwf_def)
        ! Ice sheet -> OBM: freshwater flux from the domain's ice mass balance
        ! (calc_fwf, on the Yelmo grid, restricted by the hydrographic mask) into
        ! the box model's northern (fn) or southern (fs) input flux.
        ! dom is intent(inout) because calc_fwf takes its mass-balance fields as
        ! non-intent (modifiable) allocatable dummies.
        type(ice_domain),      intent(inout) :: dom
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

    subroutine coupling_obm2ism(dom, obm, obm_name, hemisphere)
        ! OBM -> ice sheet: broadcast the box-model ocean temperature (northern box
        ! tn / southern box ts) into the domain's snapclim ocean-temperature field
        ! to_ann, which marine_shelf then reads as its ocean forcing.
        type(ice_domain), intent(inout) :: dom
        type(obm_class),  intent(in)    :: obm
        character(len=*), intent(in)    :: obm_name
        character(len=*), intent(in)    :: hemisphere

        character(len=512) :: name_loc

        name_loc = obm_name
        select case(trim(hemisphere))
            case("north")
                call calc_ocean_temperature_field(dom%snp%now%to_ann, obm%tn, name_loc)
            case("south")
                call calc_ocean_temperature_field(dom%snp%now%to_ann, obm%ts, name_loc)
        end select
    end subroutine coupling_obm2ism

    subroutine update_bipolar_hyster_forcing(t, t0, obm, dt, branch_time_thr, &
                                             rate, forcing, forc_method)
        ! Hysteresis forcing for the nautilus box model: nudge one obm control
        ! (phit/phin/fs/fn, or fn+fs) along a prescribed path. Ported verbatim
        ! from yelmox_bipolar.f90.
        real(wp),          intent(in)    :: t, t0
        type(obm_class),   intent(inout) :: obm
        real(wp),          intent(in)    :: dt
        real(wp),          intent(in)    :: branch_time_thr, rate
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

    ! ----- private helpers (ported from yelmox_bipolar.f90) -----

    function series_interp(series_time, series_var, time) result(var)
        ! Linear interpolation of a (time, var) series at `time`.
        real(wp), dimension(:), intent(in) :: series_time, series_var
        real(wp),               intent(in) :: time
        real(wp) :: var
        var = interp_linear(series_time, series_var, xout=time)
    end function series_interp

    function interp_linear(x, y, xout) result(yout)
        ! Simple linear interpolation of a point, clamped to the series endpoints.
        real(wp), dimension(:), intent(in) :: x, y
        real(wp),               intent(in) :: xout
        real(wp) :: yout
        integer  :: j, n
        real(wp) :: alph

        n = size(x)
        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do
            if (j .eq. 1) then
                yout = y(1)
            else if (j .eq. n+1) then
                yout = y(n)
            else
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if
        end if
    end function interp_linear

    function r8_normal_ab(a, b) result(val)
        ! Sample of a normal PDF with mean a, standard deviation b (Box-Muller).
        ! Ported from yelmox_bipolar.f90 (John Burkardt, MIT license).
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
