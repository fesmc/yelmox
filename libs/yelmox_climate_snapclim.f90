module yelmox_climate
    ! Climate backend adapter — snapclim variant.
    !
    ! Presents a backend-agnostic interface (climate_init / climate_update) over the
    ! legacy `snapclim`, filling a `climate_out_class`. The snapesm variant
    ! (yelmox_climate_snapesm.f90) defines the SAME module name and public API;
    ! config/Makefile_yelmox.mk compiles exactly one (CLIMATE = snapclim | snapesm).
    !
    ! The backend object is exposed as `%snp` for the one caller that reaches into
    ! snapclim internals — yelmox_bipolar/obm_coupling (only built with this variant).

    use precision,   only : wp
    use climate_out, only : climate_out_class
    use snapclim,    only : snapclim_class, snapclim_init, snapclim_update

    implicit none

    private

    type yelmox_climate_class
        type(snapclim_class) :: snp
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update

contains

    subroutine climate_init(cl, filename, domain, grid_name, nx, ny, time, basins, group)
        ! `time` is accepted for a uniform adapter signature; snapclim does not use it.
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: filename, domain, grid_name
        integer,          intent(in) :: nx, ny
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in), optional :: group

        call snapclim_init(cl%snp, filename, domain, grid_name, nx, ny, basins, group=group)

        return
    end subroutine climate_init

    subroutine climate_update(cl, out, z_srf, time, domain, dTa, dTo, dSo, dx, basins)
        type(yelmox_climate_class), intent(inout) :: cl
        type(climate_out_class),    intent(inout) :: out
        real(wp),         intent(in) :: z_srf(:,:)
        real(wp),         intent(in) :: time
        character(len=*), intent(in) :: domain
        real(wp),         intent(in), optional :: dTa, dTo, dSo, dx
        real(wp),         intent(in) :: basins(:,:)

        call snapclim_update(cl%snp, z_srf=z_srf, time=time, domain=domain, &
                             dTa=dTa, dTo=dTo, dSo=dSo, dx=dx, basins=basins)

        ! Fill the agnostic output. snapclim's reference climate is clim0.
        ! (snapclim_state_class is not exported, so copy field-by-field here.)
        out%now%tas     = cl%snp%now%tas
        out%now%pr      = cl%snp%now%pr
        out%now%tsl_ann = cl%snp%now%tsl_ann
        out%now%ta_ann  = cl%snp%now%ta_ann
        out%now%pr_ann  = cl%snp%now%pr_ann
        out%now%to_ann  = cl%snp%now%to_ann
        out%now%so_ann  = cl%snp%now%so_ann
        out%now%depth   = cl%snp%now%depth

        out%ref%tas     = cl%snp%clim0%tas
        out%ref%pr      = cl%snp%clim0%pr
        out%ref%tsl_ann = cl%snp%clim0%tsl_ann
        out%ref%ta_ann  = cl%snp%clim0%ta_ann
        out%ref%pr_ann  = cl%snp%clim0%pr_ann
        out%ref%to_ann  = cl%snp%clim0%to_ann
        out%ref%so_ann  = cl%snp%clim0%so_ann
        out%ref%depth   = cl%snp%clim0%depth

        return
    end subroutine climate_update

end module yelmox_climate
