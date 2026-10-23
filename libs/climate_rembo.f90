module climate_rembo
    ! REMBO (rembo1) as the atmosphere and surface mass balance of the rembo
    ! climate backend (yelmox_climate). Built with `make ... rembo=1`; otherwise
    ! climate_rembo_stub.f90 provides this interface and stops at init.
    !
    ! REMBO keeps its state module-global (rembo_ann), reads its own parameters
    ! (options_rembo) and runs on its compiled grid, which must be the climate
    ! grid. It updates its energy balance and surface mass balance on its own
    ! intervals (dtime_emb, dtime_smb).

    use precision,      only : wp
    use rembo_sclimate, only : rembo_init, rembo_update, rembo_equilibrate, &
                               rembo_ann, rembo_restart_write, rembo_dT_ann

    implicit none

    private

    integer, parameter :: dp = kind(1.0d0)

    public :: rembo_clim_init
    public :: rembo_clim_update
    public :: rembo_clim_dT_ann
    public :: rembo_clim_restart_write

contains

    subroutine rembo_clim_init(time, nx, ny)
        ! Initialize REMBO and check that its grid is the climate grid (nx, ny).
        real(dp), intent(in) :: time
        integer,  intent(in) :: nx, ny

        call rembo_init(time)

        if (size(rembo_ann%smb,1) /= nx .or. size(rembo_ann%smb,2) /= ny) then
            write(*,*) "rembo_clim_init:: error: REMBO's grid must be the climate grid (grid_clim)."
            write(*,*) "  REMBO:     ", size(rembo_ann%smb,1), size(rembo_ann%smb,2)
            write(*,*) "  grid_clim: ", nx, ny
            error stop 1
        end if

    end subroutine rembo_clim_init

    subroutine rembo_clim_update(time, time_ins, dT_summer, z_srf, H_ice, z_sl, equil, &
                                 ta_ann, ta_sum, pr_ann, smb, tsrf)
        ! Advance REMBO with the geometry on its grid, or, with equil, equilibrate
        ! it (cold start); return its annual and summer near-surface air
        ! temperature [K], precipitation [mm/a], surface mass balance [mm/a w.e.]
        ! and surface temperature [K].
        real(dp), intent(in)  :: time            ! [yr] model time
        real(dp), intent(in)  :: time_ins        ! [yr] time of the insolation
        real(wp), intent(in)  :: dT_summer       ! [K] summer temperature anomaly
        real(wp), intent(in)  :: z_srf(:,:), H_ice(:,:), z_sl(:,:)
        logical,  intent(in)  :: equil
        real(wp), intent(out) :: ta_ann(:,:), ta_sum(:,:), pr_ann(:,:), smb(:,:), tsrf(:,:)

        if (equil) then
            call rembo_equilibrate(time, real(z_srf, dp), real(H_ice, dp), &
                                   real(z_sl, dp), time_tot=10.0_dp)
        else
            call rembo_update(time, time_ins, real(dT_summer, dp), &
                              real(z_srf, dp), real(H_ice, dp), real(z_sl, dp))
        end if

        ta_ann = real(rembo_ann%T_ann, wp)
        ta_sum = real(rembo_ann%T_jja, wp)
        pr_ann = real(rembo_ann%pr,    wp)
        smb    = real(rembo_ann%smb,   wp)
        tsrf   = real(rembo_ann%T_srf, wp)

    end subroutine rembo_clim_update

    function rembo_clim_dT_ann(dT_summer) result(dT_ann)
        ! The annual mean air temperature anomaly [K] REMBO applies for the
        ! summer anomaly dT_summer [K] (its seasonal cycle, T_wintfac of its
        ! parameters). After rembo_clim_init.
        real(wp), intent(in) :: dT_summer
        real(wp) :: dT_ann

        dT_ann = real(rembo_dT_ann(real(dT_summer, dp)), wp)

    end function rembo_clim_dT_ann

    subroutine rembo_clim_restart_write(filename, time, z_srf, H_ice, z_sl)
        ! REMBO's restart, with the geometry on its grid.
        character(len=*), intent(in) :: filename
        real(dp),         intent(in) :: time
        real(wp),         intent(in) :: z_srf(:,:), H_ice(:,:), z_sl(:,:)

        call rembo_restart_write(filename, time, real(z_srf, dp), &
                                 real(H_ice, dp), real(z_sl, dp))

    end subroutine rembo_clim_restart_write

end module climate_rembo
