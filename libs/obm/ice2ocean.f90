module ice2ocean

    implicit none

    private
    public :: calc_fwf

    integer,  parameter :: wp = kind(1.0)


contains

function calc_fwf(rho_water,rho_ice,sec_year,mb,smb,bmb,cmb,Hice,dHidt,f_grnd,dx,dy,mask,hemisphere, fwf_def) result(fwf)
    ! Freshwater flux [Sv] from an ice sheet into the ocean box model, from
    ! the Yelmo mass balance (applied fields, m ice/yr) summed over the
    ! hydrographic mask (mask == 1: cells draining into the box), except dVdt.
    ! Positive for a mass loss of the ice sheet. Definitions (fwf_def):
    !   "dVdt"     : -dV/dt, the change in ice volume of the whole domain (dHidt)
    !   "mb"       : -net mass balance applied (mb_net)
    !   "ablation" : -the negative parts of smb, bmb and cmb (ablation only)
    ! Every definition converts ice volume to water volume (rho_ice/rho_water).
    ! The inputs are not modified. Hice, f_grnd: unused (kept for the
    ! interface).

    implicit none

    real(wp),         intent(in) :: rho_water, rho_ice   ! [kg m-3]
    real(wp),         intent(in) :: sec_year             ! [s yr-1]
    real(wp),         intent(in) :: mb(:,:), smb(:,:), bmb(:,:), cmb(:,:)   ! [m ice/yr]
    real(wp),         intent(in) :: Hice(:,:)            ! [m]
    real(wp),         intent(in) :: dHidt(:,:)           ! [m ice/yr]
    real(wp),         intent(in) :: f_grnd(:,:)          ! [1]
    real(wp),         intent(in) :: dx, dy               ! [m]
    real(wp),         intent(in) :: mask(:,:)            ! hydrographic mask (1: in the basin)
    character(len=*), intent(in) :: hemisphere
    character(len=*), intent(in) :: fwf_def
    real(wp) :: fwf                                      ! [Sv]

    ! Local variables
    logical, allocatable :: mask_basin(:,:)
    real(wp) :: dV_ice                                   ! [m3 ice/yr]
    real(wp) :: conv_m3a_Sv

    ! Conversion: ice volume [m3/yr] -> freshwater flux [Sv]. Ice volume ->
    ! the equivalent water volume is rho_ice/rho_water; /sec_year gives m3/s
    ! and 1e-6 gives Sv. Matches yelmo_regions.f90.
    conv_m3a_Sv = 1e-6_wp*(rho_ice/rho_water)/sec_year

    mask_basin = (mask .eq. 1.0_wp)

    ! Ice volume change rate [m3 ice/yr] (negative: loss)
    select case(trim(fwf_def))
        case("dVdt")
            dV_ice = sum(dHidt)*dx*dy
        case("mb")
            dV_ice = sum(mb, mask=mask_basin)*dx*dy
        case("ablation")
            dV_ice = ( sum(min(smb,0.0_wp), mask=mask_basin) &
                     + sum(min(bmb,0.0_wp), mask=mask_basin) &
                     + sum(min(cmb,0.0_wp), mask=mask_basin) )*dx*dy
        case DEFAULT
            write(*,*) "calc_fwf:: error: fwf_definition not recognized: ", trim(fwf_def)
            write(*,*) "  hemisphere: ", trim(hemisphere)
            write(*,*) "  Options: dVdt, mb, ablation."
            error stop 1
    end select

    fwf = -dV_ice*conv_m3a_Sv

    return

end function calc_fwf

end module ice2ocean
