module climate_out
    ! Climate-model-agnostic boundary-forcing output. A backend adapter
    ! (yelmox_climate) fills `now` and `ref` after each update; the domain and its
    ! coupling read only this struct, so nothing downstream depends on which
    ! backend produced the forcing. Holding both `now` and `ref` lets output
    ! writers form anomalies (now − ref) directly.
    !
    ! Every backend supplies the atmosphere. The ocean comes either as depth
    ! profiles (to_ann/so_ann, interpolated by the marine shelf) or already at the
    ! shelf base (T_shlf/S_shlf, has_ocn_shelf). A backend may also supply the
    ! surface mass balance directly (has_smb) and subglacial discharge (has_Qd).

    use precision, only : wp

    implicit none

    private

    type clim_state_class
        ! atmosphere
        real(wp), allocatable :: tas(:,:,:)     ! near-surface air temperature, monthly [K]
        real(wp), allocatable :: pr(:,:,:)      ! precipitation, monthly [mm/d]
        real(wp), allocatable :: tsl_ann(:,:)   ! sea-level air temperature, annual [K]
        real(wp), allocatable :: ta_ann(:,:)    ! near-surface air temperature, annual [K]
        real(wp), allocatable :: pr_ann(:,:)    ! precipitation, annual [mm/a]
        ! ocean as depth profiles
        real(wp), allocatable :: to_ann(:,:,:)  ! ocean temperature over depth, annual [K]
        real(wp), allocatable :: so_ann(:,:,:)  ! ocean salinity over depth, annual [psu]
        real(wp), allocatable :: depth(:)       ! ocean depth axis [m]
        ! ocean at the shelf base
        real(wp), allocatable :: T_shlf(:,:)    ! ocean temperature [K]
        real(wp), allocatable :: S_shlf(:,:)    ! ocean salinity [psu]
        real(wp), allocatable :: dT_shlf(:,:)   ! temperature anomaly to the reference ocean [K]
        real(wp), allocatable :: dS_shlf(:,:)   ! salinity anomaly to the reference ocean [psu]
        ! surface mass balance, at the current surface
        real(wp), allocatable :: smb(:,:)       ! annual (units of smbpal's smb)
        real(wp), allocatable :: tsrf(:,:)      ! surface temperature, annual [K]
        ! subglacial discharge
        real(wp), allocatable :: Qd(:,:)        ! annual (units of Yelmo's bnd%Qd)
    end type clim_state_class

    type climate_out_class
        type(clim_state_class) :: now           ! current climate
        type(clim_state_class) :: ref           ! reference climate (anomaly baseline)
        logical :: has_ocn_shelf = .false.      ! ocean given at the shelf base (now%T_shlf...)
        logical :: has_smb       = .false.      ! surface mass balance given (now%smb...)
        logical :: has_Qd        = .false.      ! subglacial discharge given (now%Qd)
    end type climate_out_class

    public :: clim_state_class
    public :: climate_out_class

end module climate_out
