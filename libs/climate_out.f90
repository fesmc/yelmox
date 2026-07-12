module climate_out
    ! Climate-model-agnostic boundary-forcing output. A backend adapter
    ! (yelmox_climate) fills `now` and `ref` after each update; the domain and its
    ! coupling read only this struct, so nothing downstream depends on whether the
    ! forcing came from snapclim or snapesm. Holding both `now` and `ref` lets output
    ! writers form anomalies (now − ref) directly.

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
        ! ocean
        real(wp), allocatable :: to_ann(:,:,:)  ! ocean temperature over depth, annual [K]
        real(wp), allocatable :: so_ann(:,:,:)  ! ocean salinity over depth, annual [psu]
        real(wp), allocatable :: depth(:)       ! ocean depth axis [m]
    end type clim_state_class

    type climate_out_class
        type(clim_state_class) :: now           ! current climate
        type(clim_state_class) :: ref           ! reference climate (anomaly baseline)
    end type climate_out_class

    public :: clim_state_class
    public :: climate_out_class

end module climate_out
