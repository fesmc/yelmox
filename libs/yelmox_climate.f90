module yelmox_climate
    ! Climate backend of a domain, chosen at runtime ([coupling] climate):
    !   "snapclim" -- snapshot/anomaly climate (snapclim)
    !   "snapesm"  -- snapshots blended by indices (snapesm)
    !
    ! Presents one interface (climate_init / climate_update) over the backends,
    ! filling a `climate_out_class`; the domain and its coupling read only that.
    ! The driver's transient forcing (tsforcing: f_now and the gains f_ta, f_to,
    ! f_so) is handed to the backend, which applies it in its own way.

    use precision,     only : wp
    use climate_out,   only : climate_out_class
    use snapclim,      only : snapclim_class, snapclim_init, snapclim_update, snapclim_air_anom
    use snapesm,       only : snapesm_class, snapesm_init, snapesm_update
    use kryos_forcing, only : tsforcing_class

    implicit none

    private

    type yelmox_climate_class
        character(len=16)    :: method = ""   ! snapclim | snapesm
        type(snapclim_class) :: snapclim
        type(snapesm_class)  :: snapesm
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update
    public :: climate_air_anom
    public :: climate_ocean_const

contains

    subroutine climate_init(cl, method, filename, domain, grid_name, nx, ny, time, basins, group)
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: method
        character(len=*), intent(in) :: filename, domain, grid_name
        integer,          intent(in) :: nx, ny
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in), optional :: group

        cl%method = trim(method)

        select case(trim(cl%method))
            case("snapclim")
                call snapclim_init(cl%snapclim, filename, domain, grid_name, nx, ny, basins, group=group)
            case("snapesm")
                call snapesm_init(cl%snapesm, filename, domain, grid_name, nx, ny, time, basins, group=group)
            case default
                write(*,*) "climate_init:: error: climate must be snapclim or snapesm; got ", trim(cl%method)
                error stop 1
        end select

    end subroutine climate_init

    subroutine climate_update(cl, out, z_srf, time, domain, dx, basins, tsf)
        ! Update the backend and fill `out`. With an active transient forcing,
        ! its spatially homogeneous anomalies (dTa = f_now*f_ta, dTo = f_now*f_to,
        ! dSo = f_now*f_so) go to the backend: snapclim uses them in place of its
        ! own index in its "anom" modes; snapesm adds them on top in every mode.
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        real(wp),         intent(in) :: z_srf(:,:)
        real(wp),         intent(in) :: time
        character(len=*), intent(in) :: domain
        real(wp),         intent(in) :: dx
        real(wp),         intent(in) :: basins(:,:)
        type(tsforcing_class), intent(in), optional :: tsf

        logical :: forced

        forced = .false.
        if (present(tsf)) forced = tsf%active

        select case(trim(cl%method))
            case("snapclim")
                if (forced) then
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, domain=domain, &
                                         dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapclim_update(cl%snapclim, z_srf=z_srf, time=time, domain=domain, &
                                         dx=dx, basins=basins)
                end if

                ! snapclim's reference climate is clim0.
                out%now%tas     = cl%snapclim%now%tas
                out%now%pr      = cl%snapclim%now%pr
                out%now%tsl_ann = cl%snapclim%now%tsl_ann
                out%now%ta_ann  = cl%snapclim%now%ta_ann
                out%now%pr_ann  = cl%snapclim%now%pr_ann
                out%now%to_ann  = cl%snapclim%now%to_ann
                out%now%so_ann  = cl%snapclim%now%so_ann
                out%now%depth   = cl%snapclim%now%depth

                out%ref%tas     = cl%snapclim%clim0%tas
                out%ref%pr      = cl%snapclim%clim0%pr
                out%ref%tsl_ann = cl%snapclim%clim0%tsl_ann
                out%ref%ta_ann  = cl%snapclim%clim0%ta_ann
                out%ref%pr_ann  = cl%snapclim%clim0%pr_ann
                out%ref%to_ann  = cl%snapclim%clim0%to_ann
                out%ref%so_ann  = cl%snapclim%clim0%so_ann
                out%ref%depth   = cl%snapclim%clim0%depth

            case("snapesm")
                if (forced) then
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, domain=domain, &
                                        dTa=tsf%dTa, dTo=tsf%dTo, dSo=tsf%dSo, dx=dx, basins=basins)
                else
                    call snapesm_update(cl%snapesm, z_srf=z_srf, time=time, domain=domain, &
                                        dx=dx, basins=basins)
                end if

                out%now%tas     = cl%snapesm%now%tas
                out%now%pr      = cl%snapesm%now%pr
                out%now%tsl_ann = cl%snapesm%now%tsl_ann
                out%now%ta_ann  = cl%snapesm%now%ta_ann
                out%now%pr_ann  = cl%snapesm%now%pr_ann
                out%now%to_ann  = cl%snapesm%now%to_ann
                out%now%so_ann  = cl%snapesm%now%so_ann
                out%now%depth   = cl%snapesm%now%depth

                out%ref%tas     = cl%snapesm%ref%tas
                out%ref%pr      = cl%snapesm%ref%pr
                out%ref%tsl_ann = cl%snapesm%ref%tsl_ann
                out%ref%ta_ann  = cl%snapesm%ref%ta_ann
                out%ref%pr_ann  = cl%snapesm%ref%pr_ann
                out%ref%to_ann  = cl%snapesm%ref%to_ann
                out%ref%so_ann  = cl%snapesm%ref%so_ann
                out%ref%depth   = cl%snapesm%ref%depth
        end select

    end subroutine climate_update

    function climate_air_anom(cl, time) result(dT)
        ! [K] The backend's air-temperature anomaly index at `time`, the forcing of
        ! the bipolar ocean box model: snapclim's index series `at` scaled by
        ! dTa_const. snapesm has no such index (its indices weight snapshots).
        type(yelmox_climate_class), intent(in) :: cl
        real(wp),                   intent(in) :: time
        real(wp) :: dT

        select case(trim(cl%method))
            case("snapclim")
                dT = snapclim_air_anom(cl%snapclim, time)
            case default
                write(*,*) "climate_air_anom:: error: climate = ", trim(cl%method), &
                           " has no air-temperature anomaly index (needed by the ocean box model)."
                error stop 1
        end select

    end function climate_air_anom

    logical function climate_ocean_const(cl) result(is_const)
        ! Is the backend's ocean held at its reference state (snapclim ocn_type
        ! = "const")? REMBO then adds its own ocean anomaly.
        type(yelmox_climate_class), intent(in) :: cl

        is_const = .false.
        if (trim(cl%method) == "snapclim") is_const = (trim(cl%snapclim%par%ocn_type) == "const")

    end function climate_ocean_const

end module yelmox_climate
