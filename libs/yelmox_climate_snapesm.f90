module yelmox_climate
    ! Climate backend adapter — snapesm variant.
    !
    ! Presents a backend-agnostic interface (climate_init / climate_update) over
    ! `snapesm`, filling a `climate_out_class`. The snapclim variant
    ! (yelmox_climate_snapclim.f90) defines the SAME module name and public API;
    ! config/Makefile_yelmox.mk compiles exactly one (CLIMATE = snapclim | snapesm).

    use precision,   only : wp
    use climate_out, only : climate_out_class
    use snapesm,     only : snapesm_class, snapesm_init, snapesm_update

    implicit none

    private

    type yelmox_climate_class
        type(snapesm_class) :: snp
    end type yelmox_climate_class

    public :: yelmox_climate_class
    public :: climate_init
    public :: climate_update

contains

    subroutine climate_init(cl, filename, domain, grid_name, nx, ny, time, basins, group)
        type(yelmox_climate_class), intent(inout) :: cl
        character(len=*), intent(in) :: filename, domain, grid_name
        integer,          intent(in) :: nx, ny
        real(wp),         intent(in) :: time
        real(wp),         intent(in) :: basins(:,:)
        character(len=*), intent(in), optional :: group

        call snapesm_init(cl%snp, filename, domain, grid_name, nx, ny, time, basins, group=group)

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

        call snapesm_update(cl%snp, z_srf=z_srf, time=time, domain=domain, &
                            dTa=dTa, dTo=dTo, dSo=dSo, dx=dx, basins=basins)

        ! Fill the agnostic output (climate_state_class is not exported).
        out%now%tas     = cl%snp%now%tas
        out%now%pr      = cl%snp%now%pr
        out%now%tsl_ann = cl%snp%now%tsl_ann
        out%now%ta_ann  = cl%snp%now%ta_ann
        out%now%pr_ann  = cl%snp%now%pr_ann
        out%now%to_ann  = cl%snp%now%to_ann
        out%now%so_ann  = cl%snp%now%so_ann
        out%now%depth   = cl%snp%now%depth

        out%ref%tas     = cl%snp%ref%tas
        out%ref%pr      = cl%snp%ref%pr
        out%ref%tsl_ann = cl%snp%ref%tsl_ann
        out%ref%ta_ann  = cl%snp%ref%ta_ann
        out%ref%pr_ann  = cl%snp%ref%pr_ann
        out%ref%to_ann  = cl%snp%ref%to_ann
        out%ref%so_ann  = cl%snp%ref%so_ann
        out%ref%depth   = cl%snp%ref%depth

        return
    end subroutine climate_update

end module yelmox_climate
