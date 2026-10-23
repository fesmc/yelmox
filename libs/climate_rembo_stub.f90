module climate_rembo
    ! Stand-in for the REMBO adapter (climate_rembo.f90) when yelmox is built
    ! without REMBO (rembo=0, the default): the same interface, stopping at
    ! init if climate = "rembo" is chosen.

    use precision, only : wp, dp

    implicit none

    private

    public :: rembo_clim_init
    public :: rembo_clim_update
    public :: rembo_clim_dT_ann
    public :: rembo_clim_restart_write

contains

    subroutine rembo_clim_init(time, nx, ny)
        real(dp), intent(in) :: time
        integer,  intent(in) :: nx, ny

        call not_built()

    end subroutine rembo_clim_init

    subroutine rembo_clim_update(time, time_ins, dT_summer, z_srf, H_ice, z_sl, equil, &
                                 ta_ann, ta_sum, pr_ann, smb, tsrf)
        real(dp), intent(in)  :: time, time_ins
        real(wp), intent(in)  :: dT_summer
        real(wp), intent(in)  :: z_srf(:,:), H_ice(:,:), z_sl(:,:)
        logical,  intent(in)  :: equil
        real(wp), intent(out) :: ta_ann(:,:), ta_sum(:,:), pr_ann(:,:), smb(:,:), tsrf(:,:)

        call not_built()

    end subroutine rembo_clim_update

    function rembo_clim_dT_ann(dT_summer) result(dT_ann)
        real(wp), intent(in) :: dT_summer
        real(wp) :: dT_ann

        dT_ann = 0.0_wp
        call not_built()

    end function rembo_clim_dT_ann

    subroutine rembo_clim_restart_write(filename, time, z_srf, H_ice, z_sl)
        character(len=*), intent(in) :: filename
        real(dp),         intent(in) :: time
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)

        call not_built()

    end subroutine rembo_clim_restart_write

    subroutine not_built()
        write(*,*) "climate_rembo:: error: climate = rembo needs yelmox built with REMBO &
                   &(make yelmox rembo=1)."
        error stop 1
    end subroutine not_built

end module climate_rembo
