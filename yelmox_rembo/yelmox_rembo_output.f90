module yelmox_rembo_output
    ! REMBO-flavor diagnostic output: the REMBO-specific fields that the shared
    ! per-module yelmox output does NOT carry. Everything generic (Yelmo 2D/1D,
    ! isostasy, marine shelf, hi-res hub) is written by the shared domain_write_*
    ! routines into yelmo.nc / yelmo_sm.nc / yelmo_ts.nc / isos.nc / mshlf.nc /
    ! htopo.nc, exactly as in the other flavors. This module adds two files:
    !   rembo.nc     -- 2D REMBO climate (tm_2D cadence)
    !   rembo_ts.nc  -- 1D hysteresis/forcing + REMBO integrated metrics (tm_1D)

    use ncio
    use yelmo
    use rembo_sclimate, only : rembo_class
    use hyster,         only : hyster_class

    implicit none
    private

    public :: rembo_write_2D_init, rembo_write_2D_step
    public :: rembo_write_1D_init, rembo_write_1D_step

contains

    subroutine rembo_write_2D_init(ylmo, filename, time_init, units)
        ! Create rembo.nc (xc/yc/time) for the 2D REMBO climate fields.
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename, units
        real(wp),          intent(IN) :: time_init

        call nc_create(filename)
        call nc_write_dim(filename,"xc",  x=real(ylmo%grd%G%x,wp)*1e-3, units="kilometers")
        call nc_write_dim(filename,"yc",  x=real(ylmo%grd%G%y,wp)*1e-3, units="kilometers")
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)
    end subroutine rembo_write_2D_init

    subroutine rembo_write_2D_step(rembo, filename, time)
        ! Append one record of the REMBO climate fields to rembo.nc.
        type(rembo_class), intent(IN) :: rembo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        integer :: ncid, n

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_time_index(filename,"time",time,ncid)
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"Ta_ann",rembo%T_ann,units="K", &
                      long_name="REMBO near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",rembo%T_jja,units="K", &
                      long_name="REMBO near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",rembo%pr*1e-3,units="m/a water equiv.", &
                      long_name="REMBO precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ann",rembo%smb*1e-3,units="m/a water equiv.", &
                      long_name="REMBO surface mass balance (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_close(ncid)
    end subroutine rembo_write_2D_step

    subroutine rembo_write_1D_init(filename, time_init, units, hyst)
        ! Create rembo_ts.nc (time + dT_axis) for the hysteresis diagnostics. The
        ! dT_axis dimension carries the ice volume as a function of forcing (V_dT).
        character(len=*),   intent(IN) :: filename, units
        real(wp),           intent(IN) :: time_init
        type(hyster_class), intent(IN) :: hyst

        integer :: n
        real(wp), allocatable :: dT_axis(:)

        call nc_create(filename)
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        allocate(dT_axis(1000))
        do n = 1, 1000
            dT_axis(n) = hyst%par%f_min + (hyst%par%f_max-hyst%par%f_min)*(n-1)/real(1000-1,wp)
        end do
        call nc_write_dim(filename,"dT_axis",x=dT_axis,units="degC")

        ! Initialize V_dT with missing values (filled per forcing bin at runtime).
        dT_axis = missing_value
        call nc_write(filename,"V_dT",dT_axis,dim1="dT_axis",missing_value=missing_value)
    end subroutine rembo_write_1D_init

    subroutine rembo_write_1D_step(ylmo, hyst, rembo, filename, time, dT_ann, dT_ocn)
        ! Append one record of the hysteresis/forcing scalars, the ice volume in
        ! forcing phase space (V_dT), and the REMBO integrated metrics (smb_mean,
        ! aar) to rembo_ts.nc. Generic 1D aggregates (V_ice, A_ice, ...) live in
        ! the shared yelmo_ts.nc.
        type(yelmo_class),  intent(IN) :: ylmo
        type(hyster_class), intent(IN) :: hyst
        type(rembo_class),  intent(IN) :: rembo
        character(len=*),   intent(IN) :: filename
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: dT_ann
        real(wp), intent(IN) :: dT_ocn

        integer  :: ncid, n, k
        real(wp) :: npmb, ntot, aar, smb_tot
        real(wp) :: dT_axis(1000)

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_time_index(filename,"time",time,ncid)
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! ===== Hysteresis / forcing scalars =====
        call nc_write(filename,"hyst_f_now",hyst%f_now,units="K",long_name="hyst: forcing value", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 yr)", &
                      long_name="hyst: forcing rate of change",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt_ave,units="m/yr", &
                      long_name="hyst: rms thickness rate of change",dim1="time",start=[n],ncid=ncid)

        ! Ice volume in forcing (dT) phase space.
        call nc_read(filename,"dT_axis",dT_axis)
        k = minloc(abs(dT_axis-hyst%f_now),dim=1)
        call nc_write(filename,"V_dT",ylmo%reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="dT_axis",start=[k],ncid=ncid)

        ! Yelmo model metrics (model speed, dt, eta).
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! ===== REMBO integrated metrics (smb_tot [Gt/yr], aar [-]) =====
        ntot = count(ylmo%tpo%now%H_ice .gt. 0.0)
        if (ntot .gt. 0.0) then
            npmb = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. rembo%smb .gt. 0.0)
            aar  = real(npmb,prec) / real(ntot,prec)
            smb_tot = (ylmo%tpo%par%dx**2)*sum(rembo%smb*1e-3,mask=ylmo%tpo%now%H_ice .gt. 0.0)
            ! [m^3/yr] * [1000 kg/m^3] * [1e-12 Gt/kg] == [Gt/yr]
            smb_tot = smb_tot * 1000 * 1e-12
        else
            aar     = 0.0
            smb_tot = 0.0
        end if

        call nc_write(filename,"dT_jja",hyst%f_now,units="K",long_name="Temp. anomaly, regional JJA mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dT_ann",dT_ann,units="K",long_name="Temp. anomaly, regional annual mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dT_ocn",dT_ocn,units="K",long_name="Temp. anomaly, regional oceanic mean", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb_mean",smb_tot,units="Gt/yr",long_name="Mean smb over the ice sheet", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"aar",aar,units="1",long_name="Accumulation area ratio", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_close(ncid)
    end subroutine rembo_write_1D_step

end module yelmox_rembo_output
