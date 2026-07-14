module yelmox_esm_output
    ! Output writers for the multigrid ESM driver (yelmox_esm).
    !
    ! Ported verbatim from yelmox_esm.f90's contained writers so the ESM driver
    ! keeps identical NetCDF output. The two 2D writers that read the (program-
    ! local) ctrl_params ESM/SMB switch now take an explicit use_smb logical
    ! instead of host-associating it.
    !
    !   write_step_2D_combined   standard heavy 2D output
    !   write_step_2D_small      small 2D output
    !   write_1D_esm             ESM 1D timeseries
    !   write_step_2D_cmip       CMIP/ISMIP7-formatted 2D output
    !   write_step_1D_cmip       CMIP/ISMIP7-formatted 1D output

    use nml
    use ncio
    use yelmo
    use esm_forcing
    use fastisostasy    ! isos_class (reexports barysealevel)
    use marine_shelf
    use smbpal

    implicit none

    private
    public :: write_step_2D_combined
    public :: write_step_2D_small
    public :: write_1D_esm
    public :: write_step_2D_cmip
    public :: write_step_1D_cmip

contains

    subroutine write_step_2D_combined(ylmo,isos,esm,mshlf,srf,use_smb,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: srf
        logical,                 intent(IN) :: use_smb

        character(len=*),       intent(IN) :: filename
        real(wp),               intent(IN) :: time

        ! Local variables
        integer  :: ncid, n

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_grz",ylmo%tpo%now%mask_grz,units="",long_name="Grounding-zone mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lsf",ylmo%tpo%now%lsf,units="",long_name="LSF mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_ice",ylmo%bnd%mask_ice,units="",long_name="Ice mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb_flt+ylmo%tpo%now%cmb_grnd,units="m/a ice equiv.",long_name="Calving mass balance rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a ice equiv.",long_name="Calving mass balance rate flt", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_grnd",ylmo%tpo%now%cmb_grnd,units="m/a ice equiv.",long_name="Calving mass balance rate grnd", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_bmb",ylmo%tpo%now%f_grnd_bmb,units="1",long_name="Grounded fraction (bmb)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
                      long_name="Distance to nearest grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity (z)", &
                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"hyd_W_til",ylmo%hyd%now%W_til,units="m",long_name="Basal water layer thickness", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/a ice equiv.",long_name="Net surface mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Net basal mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Net margin-front mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! External data
        call nc_write(filename,"dzbdt",isos%out%dwdt,units="m/a",long_name="Bedrock uplift rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
    
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! ESM Atmospheric boundary fields            
        call nc_write(filename,"t2m_ann",esm%t2m_ann+SUM(esm%dts, dim=3)/12.0,units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"t2m_sum",esm%t2m_sum+0.333*(esm%dts(:,:,12)+esm%dts(:,:,1)+esm%dts(:,:,2)),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_ann",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        if (use_smb) then
            call nc_write(filename,"dsmb_ann",1e-3*SUM(esm%dsmb, dim=3)/12.0,units="m/a water equiv.",long_name="SMB anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else
            call nc_write(filename,"pr_ann",esm%pr_ann*1e-3*esm%dpr(:,:,1),units="m/a water equiv.",long_name="Precipitation (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            !call nc_write(filename,"dpr_var",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Precipitation anomaly (variability)", &
            !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        end if
        
        ! Oceanic boundary conditions
        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto",esm%dto,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso",esm%dso,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"dto_var",esm%dto_var,units="K",long_name="Shelf temperature anomaly (variability)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"dso_var",esm%dso_var,units="PSU",long_name="Shelf salinity anomaly (variability)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dS_shlf",mshlf%now%dS_shlf,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_fp_shlf",mshlf%now%T_fp_shlf,units="K",long_name="Shelf freezing temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="", &
                        long_name="Ocean mask (0: land, 1: grline, 2: fltline, 3: open ocean, 4: deep ocean, 5: lakes)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"tf_basin",mshlf%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr_basin",mshlf%now%tf_corr_basin,units="K",long_name="Shelf thermal forcing basin-wide correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        if (trim(mshlf%par%bmb_method) .eq. "pico") then
            call nc_write(filename,"d_shlf",mshlf%pico%now%d_shlf,units="km",long_name="Shelf distance to grounding line", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"d_if",mshlf%pico%now%d_if,units="km",long_name="Shelf distance to ice front", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"boxes",mshlf%pico%now%boxes,units="",long_name="Shelf boxes", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"r_shlf",mshlf%pico%now%r_shlf,units="",long_name="Ratio of ice shelf", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"T_box",mshlf%pico%now%T_box,units="K?",long_name="Temperature of boxes", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"S_box",mshlf%pico%now%S_box,units="PSU",long_name="Salinity of boxes", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"A_box",mshlf%pico%now%A_box*1e-6,units="km2",long_name="Box area of ice shelf", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        end if

        call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! Strain-rate and stress tensors 
        if (.FALSE.) then

            call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
                          dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
                          dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"eps_eig_1",ylmo%mat%now%strn2D%eps_eig_1,units="1/yr",long_name="Eigen strain 1", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eig_2",ylmo%mat%now%strn2D%eps_eig_2,units="1/yr",long_name="Eigen strain 2", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eff",ylmo%tpo%now%eps_eff,units="yr^-1",long_name="Effective calving strain", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"tau_eig_1",ylmo%mat%now%strs2D%tau_eig_1,units="Pa",long_name="Eigen stress 1", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eig_2",ylmo%mat%now%strs2D%tau_eig_2,units="Pa",long_name="Eigen stress 2", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eff",ylmo%tpo%now%tau_eff,units="Pa",long_name="Effective calving stress", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_combined

    subroutine write_step_2D_small(ylmo,isos,esm,mshlf,smbp,use_smb,filename,time)

        implicit none

        type(yelmo_class),       intent(IN) :: ylmo
        type(isos_class),        intent(IN) :: isos
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        type(smbpal_class),      intent(IN) :: smbp
        logical,                 intent(IN) :: use_smb

        character(len=*),        intent(IN) :: filename
        real(wp),                intent(IN) :: time

        ! Local variables
        integer  :: ncid, n

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)    
        call nc_write(filename,"cb_tgt",ylmo%dyn%now%cb_tgt,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! === yelmo forcing ===
        ! ESM Atmospheric boundary fields            
        call nc_write(filename,"t2m_ann",esm%t2m_ann+esm%dts(:,:,1),units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"t2m_sum",esm%t2m_sum+esm%dts(:,:,1),units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts",SUM(esm%dts, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dts_var",SUM(esm%dts_var, dim=3)/12.0,units="K",long_name="Surface air temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ann",ylmo%tpo%now%smb,units="m/a water equiv.",long_name="SMB (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        if (use_smb) then
            call nc_write(filename,"dsmb_ann",1e-3*SUM(esm%dsmb, dim=3)/12.0,units="m/a water equiv.",long_name="SMB anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dsmbdz",1e-3*esm%dsmbdz,units="m/a m-1 water equiv.",long_name="SMB lapse rate", &
                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else
            call nc_write(filename,"pr_ann",esm%pr_ann*1e-3*esm%dpr(:,:,1),units="m/a water equiv.",long_name="Precipitation (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_ann",SUM(esm%dpr, dim=3)/12.0,units="%",long_name="Precipitation anomaly (ann)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"dpr_var",SUM(esm%dpr_var, dim=3)/12.0,units="%",long_name="Precipitation anomaly (variability)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        end if

        ! Oceanic boundary conditions
        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto",esm%dto,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso",esm%dso,units="PSU",long_name="Shelf salinity anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dto_var",esm%dto_var,units="K",long_name="Shelf temperature anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dso_var",esm%dso_var,units="PSU",long_name="Shelf salinity anomaly (variability)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        if (.FALSE.) then
            call nc_write(filename,"so_ref",esm%so_ref%var(:,:,1,1),units="PSU",long_name="Reference oceanic salinity", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"to_ref",esm%to_ref%var(:,:,1,1),units="K",long_name="Reference oceanic temperature", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)    
            call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/a ice equiv.",long_name="Net surface mass balance", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)    
        end if
        call nc_write(filename,"Qd_ann",esm%Qd_ann,units="m3/s",long_name="Subglacial discharge (annual)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Qd_sum",esm%Qd_sum,units="m3/s",long_name="Sunglacial discharge (summer)", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_small

    ! ===== esm output routines =========
    subroutine write_1D_esm(ylmo, esm, mshlf, filename, time)

        ! Used to plot climatic variable fields
    
        implicit none
    
        type(yelmo_class),       intent(IN) :: ylmo
        type(esm_forcing_class), intent(IN) :: esm
        type(marshelf_class),    intent(IN) :: mshlf
        character(len=*),        intent(IN) :: filename
        real(wp),                intent(IN) :: time
    
        ! Local variables
        type(yregions_class) :: reg
    
        integer  :: ncid, n
        real(wp) :: rho_ice, density_corr, m3yr_to_kgs, esm_correction, yr_to_sec
    
        real(wp) :: dx, dy
        integer  :: npts_tot, npts_flt
        real(wp) :: smb_tot, bmb_shlf_t
    
        ! Climatic variables - atmosphere
        real(wp) :: t2m_1d, pr_1d
        real(wp) :: dt_1d,  dt_var_1d, dpr_1d, dpr_var_1d
    
        ! Climatic variables - ocean
        real(wp) :: to_1d
        real(wp) :: so_1d
        real(wp) :: tf_1d
        real(wp) :: dto_1d, dto_var_1d
        real(wp) :: dso_1d, dso_var_1d

        ! Missing-value (N/A) sentinel for the ocean shelf-draft diagnostics.
        ! The marine-shelf module leaves a large out-of-range fill value in cells
        ! where no ocean data is available (e.g. a minority of Greenland floating
        ! cells, whose tf forcing carries a NaN _FillValue). ocn_mean below drops
        ! those PER FIELD (fill is not identical across T_shlf/S_shlf/tf_shlf) via
        ! a physical bound, returning this sentinel only if a field has no valid
        ! floating-ice cell.
        real(wp), parameter :: mv_ocn = -9999.0_wp

        logical, allocatable :: mask_tot(:,:)
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:)
    
        dx = ylmo%grd%G%dx
        dy = ylmo%grd%G%dy
    
        allocate(mask_tot (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grnd(ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_flt (ylmo%grd%G%nx, ylmo%grd%G%ny))
    
        ! === Unit conversion factors =========================================
        rho_ice        = 917.0_wp           ! ice density kg m-3
        m3yr_to_kgs    = 3.2e-5_wp          ! m3 yr-1 pure water -> kg s-1
        density_corr   = rho_ice / 1000.0_wp
        esm_correction = m3yr_to_kgs * density_corr
        yr_to_sec      = 31556952.0_wp
    
        ! === Masks ===========================================================
    
        mask_tot  = (ylmo%tpo%now%H_ice .gt. 0.0_wp)
        mask_grnd = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp)
        mask_flt  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
    
        npts_tot = count(mask_tot)
        npts_flt = count(mask_flt)
    
        ! === Regional object =================================================
    
        reg = ylmo%reg
    
        ! === Integrated fluxes [m yr-1 * m2 -> m3 yr-1] =====================
        ! Total SMB over all ice-covered cells [m3 yr-1]
        smb_tot    = sum(ylmo%bnd%smb,       mask=mask_tot)  * (dx * dy)
    
        ! Total BMB beneath floating ice [m3 yr-1]
        bmb_shlf_t = sum(ylmo%bnd%bmb_shlf, mask=mask_flt)  * (dx * dy)
    
        ! === Spatially averaged climatic fields ==============================
    
        ! Atmosphere (averaged over all ice)
        if (npts_tot .gt. 0.0) then
            t2m_1d     = sum(esm%t2m_ann + esm%dts(:,:,1),      mask=mask_tot) / npts_tot
            pr_1d      = sum(esm%pr_ann * 1e-3_wp * esm%dpr(:,:,1), mask=mask_tot) / npts_tot
            dt_1d      = sum(esm%dts(:,:,1),                     mask=mask_tot) / npts_tot
            dpr_1d     = sum(100.0_wp * esm%dpr(:,:,1),          mask=mask_tot) / npts_tot
            dt_var_1d  = sum(esm%dts_var(:,:,1),                 mask=mask_tot) / npts_tot
            dpr_var_1d = sum(100.0_wp * esm%dpr_var(:,:,1),      mask=mask_tot) / npts_tot
        else
            t2m_1d = 0.0_wp; pr_1d = 1.0_wp; dt_1d = 0.0_wp; dpr_1d = 1.0_wp
            dt_var_1d = 0.0_wp; dpr_var_1d = 0.0_wp
        end if

        ! Ocean (averaged over floating ice, per field, excluding fill cells).
        ! Each shelf field is averaged only over floating cells whose value is
        ! finite and inside a physical bound, so marine-shelf fill cells (large
        ! out-of-range sentinel) no longer contaminate the mean. A field returns
        ! mv_ocn (N/A) if it has no valid floating cell -- e.g. Greenland tf_1d,
        ! whose tf_shlf (~273 K here, an artifact of feeding tf-as-temperature)
        ! falls outside the plausible thermal-forcing range; the Greenland shelf
        ! temperature/salinity (to_1d/so_1d) remain valid.
        to_1d      = ocn_mean(mshlf%now%T_shlf,  mask_flt,  240.0_wp, 320.0_wp, mv_ocn)
        so_1d      = ocn_mean(mshlf%now%S_shlf,  mask_flt,    0.0_wp,  60.0_wp, mv_ocn)
        tf_1d      = ocn_mean(mshlf%now%tf_shlf, mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dto_1d     = ocn_mean(esm%dto,           mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dso_1d     = ocn_mean(esm%dso,           mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dto_var_1d = ocn_mean(esm%dto_var,       mask_flt, -100.0_wp, 100.0_wp, mv_ocn)
        dso_var_1d = ocn_mean(esm%dso_var,       mask_flt, -100.0_wp, 100.0_wp, mv_ocn)

        ! === Write to file ===================================================
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
            
        ! -- Variability fields -----------------------------------------------
        call nc_write(filename, "dt_var_1d",  dt_var_1d,  units="K",   &
            long_name="Mean ice surf. Temp. Anomaly (Variability)",     &
            standard_name="Mean ice surf. Temp. Anomaly (Variability)", &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dpr_var_1d", dpr_var_1d, units="%",   &
            long_name="Mean ice surf. Pr. Anomaly (Variability)",       &
            standard_name="Mean ice surf. Pr. Anomaly (Variability)",   &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dto_var_1d", dto_var_1d, units="K",   &
            long_name="Mean ice-shelf Temp. Anomaly (Variability)",     &
            standard_name="Mean ice-shelf Temp. Anomaly (Variability)", &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dso_var_1d", dso_var_1d, units="PSU", &
            long_name="Mean ice-shelf Sal. Anomaly (Variability)",      &
            standard_name="Mean ice-shelf Sal. Anomaly (Variability)",  &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Atmosphere fields ------------------------------------------------
        call nc_write(filename, "t2m_1d",  t2m_1d,  units="K",       &
            long_name="Mean ice surf. Temp.",                          &
            standard_name="Mean ice surf. Temp.",                      &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "pr_1d",   pr_1d,   units="m yr-1",  &
            long_name="Mean ice surf. Pr.",                            &
            standard_name="Mean ice surf. Pr.",                        &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dt_1d",   dt_1d,   units="K",       &
            long_name="Mean ice surf. Temp. Anomaly",                  &
            standard_name="Mean ice surf. Temp. Anomaly",              &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dpr_1d",  dpr_1d,  units="%",       &
            long_name="Mean ice surf. Pr. Anomaly",                    &
            standard_name="Mean ice surf. Pr. Anomaly",                &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Ocean fields -----------------------------------------------------
        call nc_write(filename, "to_1d",   to_1d,   units="K",       &
            long_name="Mean ice-shelf Temp.",                          &
            standard_name="Mean ice-shelf draft Temp.",                &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "so_1d",   so_1d,   units="PSU",     &
            long_name="Mean ice-shelf draft Sal.",                     &
            standard_name="Mean ice-shelf draft Sal.",                 &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "tf_1d",   tf_1d,   units="K",       &
            long_name="Mean ice-shelf TF.",                            &
            standard_name="Mean ice-shelf draft TF",                   &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dto_1d",  dto_1d,  units="K",       &
            long_name="Mean ice-shelf Temp. Anomaly",                  &
            standard_name="Mean ice-shelf Temp. Anomaly",              &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "dso_1d",  dso_1d,  units="PSU",     &
            long_name="Mean ice-shelf draft Sal. Anomaly",             &
            standard_name="Mean ice-shelf draft Sal. Anomaly",         &
            dim1="time", start=[n], ncid=ncid)
    
        ! -- Integrated mass fluxes [kg s-1] ----------------------------------
        call nc_write(filename, "smb_tot",  smb_tot  * esm_correction, units="kg s-1", &
            long_name="Total SMB flux",                                                 &
            standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",     &
            dim1="time", start=[n], ncid=ncid)
        call nc_write(filename, "bmb_shlf", bmb_shlf_t * esm_correction, units="kg s-1", &
            long_name="Total BMB flux beneath floating ice",                              &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",         &
            dim1="time", start=[n], ncid=ncid)
    
        call nc_close(ncid)
    
        return
    
    end subroutine write_1D_esm

    function ocn_mean(field, base_mask, lo, hi, mv) result(val)
        ! Mean of `field` over the `base_mask` cells whose value is finite and
        ! inside the physical range (lo,hi); returns `mv` (missing value / N/A)
        ! if no cell qualifies. NaN/Inf fail the comparisons and so are dropped.
        ! Keeps marine-shelf fill cells (large out-of-range sentinel) out of the
        ! floating-ice ocean-forcing diagnostics, per field.
        implicit none
        real(wp), intent(IN) :: field(:,:)
        logical,  intent(IN) :: base_mask(:,:)
        real(wp), intent(IN) :: lo, hi, mv
        real(wp) :: val
        ! Local variables
        logical, allocatable :: m(:,:)
        integer :: npts

        allocate(m(size(field,1), size(field,2)))
        m = base_mask .and. field .gt. lo .and. field .lt. hi
        npts = count(m)
        if (npts .gt. 0) then
            val = sum(field, mask=m) / real(npts, wp)
        else
            val = mv
        end if

        return
    end function ocn_mean

    subroutine write_step_2D_cmip(ylmo, mshlf, filename, time)
        ! Writes all mandatory (and key optional) 2-D ISMIP7 variables.
        ! ST = snapshot (end-of-year); FL = yearly-average flux.
        ! -------------------------------------------------------------------------
        
        implicit none
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        character(len=*),     intent(IN) :: filename
        real(wp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        integer  :: ncid, n, i, j, k, nz
        
        real(wp) :: rho_ice
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: esm_correction   ! [kg m-2 s-1] per [m yr-1]
        real(wp) :: yr_to_sec
        
        ! 2-D working arrays
        real(wp), allocatable :: bmb_grnd_masked(:,:), bmb_shlf_masked(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: T_top_ice(:,:), T_base_grnd(:,:), T_base_flt(:,:), T_avg(:,:)
        real(wp), allocatable :: dTdz_base_grnd(:,:), dTdz_base_flt(:,:)
        real(wp), allocatable :: flux_grl_2d(:,:), flux_clv_2d(:,:), tfbase(:,:)
        real(wp), allocatable :: ux_aa(:,:), uy_aa(:,:) 
        real(wp), allocatable :: uz_s_masked(:,:), uz_b_masked(:,:)   
        
        ! ---- allocate -------------------------------------------------------
        allocate(bmb_grnd_masked (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(bmb_shlf_masked (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(z_base          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_top_ice       (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_grnd     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_flt      (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_avg           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_grnd  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_flt   (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(flux_grl_2d     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(flux_clv_2d     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(tfbase          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(ux_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uy_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_s_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_b_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        
        ! ---- initialise -----------------------------------------------------
        bmb_grnd_masked = 0.0_wp;  bmb_shlf_masked = 0.0_wp
        z_base          = 0.0_wp;  T_top_ice        = 0.0_wp
        T_base_grnd     = 0.0_wp;  T_base_flt       = 0.0_wp
        T_avg           = 0.0_wp;  tfbase           = 0.0_wp
        dTdz_base_grnd  = 0.0_wp;  dTdz_base_flt    = 0.0_wp
        flux_grl_2d     = 0.0_wp;  flux_clv_2d      = 0.0_wp  
        ux_aa           = 0.0_wp;  uy_aa            = 0.0_wp
        
        ! ---- unit conversion ------------------------------------------------
        rho_ice        = 917.0_wp
        m3yr_to_kgs    = 3.2e-5_wp
        density_corr   = rho_ice / 1000.0_wp
        esm_correction = m3yr_to_kgs * density_corr
        yr_to_sec      = 31556952.0_wp
        nz = ylmo%dyn%par%nz_aa
        
        ! ---- derived fields -------------------------------------------------
        
        ! Ice-base elevation
        z_base = ylmo%tpo%now%z_base
        
        ! BMB masked to grounded / floating
        where (ylmo%tpo%now%f_grnd .gt. 0.0_wp) bmb_grnd_masked = ylmo%thrm%now%bmb_grnd
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) bmb_shlf_masked = ylmo%bnd%bmb_shlf
        
        ! Temperature fields
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp) T_top_ice = ylmo%thrm%now%T_ice(:,:,nz)
        where (ylmo%tpo%now%f_grnd .gt. 0.0_wp) T_base_grnd = ylmo%thrm%now%T_ice(:,:,1)
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) T_base_flt = ylmo%thrm%now%T_ice(:,:,1)
        
        ! Depth-averaged temperature
        if (nz .gt. 0) then
            do k = 1, nz
                where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
                    T_avg = T_avg + ylmo%thrm%now%T_ice(:,:,k)
                end where
            end do
            where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
                T_avg = T_avg / real(nz, wp)
            end where
        end if
        
        ! Vertical basal temperature gradient (first-order upward difference)
        if (nz .gt. 1) then
            where (ylmo%tpo%now%f_grnd .gt. 0.0_wp .and. ylmo%tpo%now%H_ice .gt. 1.0_wp)
                dTdz_base_grnd = (ylmo%thrm%now%T_ice(:,:,2) - ylmo%thrm%now%T_ice(:,:,1)) &
                                 / (ylmo%tpo%now%H_ice / real(nz - 1, wp))
            end where
            where (ylmo%tpo%now%H_ice .gt. 1.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
                dTdz_base_flt  = (ylmo%thrm%now%T_ice(:,:,2) - ylmo%thrm%now%T_ice(:,:,1)) &
                                 / (ylmo%tpo%now%H_ice / real(nz - 1, wp))
            end where
        end if
        
        ! Vertical velocities (weird shape with regions)
        uz_s_masked = 0.0_wp
        uz_b_masked = 0.0_wp
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
            uz_s_masked = ylmo%dyn%now%uz_s / yr_to_sec
            uz_b_masked = ylmo%dyn%now%uz_b / yr_to_sec
        end where

        ! Grounding-line flux (2-D, for ligroundf field)
        ! Use mask_grz convention from yelmo_calving.f90: grounded + mask_grz==0
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp &
            .and. ylmo%tpo%now%mask_grz .eq. 0.0_wp)
            flux_grl_2d = ylmo%dyn%now%uxy_bar * ylmo%tpo%now%H_ice * rho_ice / yr_to_sec
        end where
        
        flux_clv_2d = (ylmo%tpo%now%cmb_flt+ylmo%tpo%now%cmb_grnd) * esm_correction

        ! Thermal forcing at ice base (floating only)
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
            tfbase = mshlf%now%tf_shlf
        end where
        
        ! Mean velocities interpolated onto aa-nodes (staggered → centred)
        do j = 2, ylmo%grd%G%ny - 1
        do i = 2, ylmo%grd%G%nx - 1
            ux_aa(i,j) = 0.5_wp * (ylmo%dyn%now%ux_bar(i,j) + ylmo%dyn%now%ux_bar(i-1,j))
            uy_aa(i,j) = 0.5_wp * (ylmo%dyn%now%uy_bar(i,j) + ylmo%dyn%now%uy_bar(i,j-1))
        end do
        end do
        
        ! ---- open file & find time index ------------------------------------
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! ====================================================================
        ! 2-D ST variables  (snapshot)
        ! ====================================================================
        
        call nc_write(filename, "lithk", ylmo%tpo%now%H_ice, &
            units="m", long_name="Ice thickness", &
            standard_name="land_ice_thickness", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "orog", ylmo%tpo%now%z_srf, &
            units="m", long_name="Surface elevation", &
            standard_name="surface_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "topg", ylmo%bnd%z_bed, &
            units="m", long_name="Bedrock elevation", &
            standard_name="bedrock_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "base", z_base, &
            units="m", long_name="Ice base elevation", &
            standard_name="base_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelsurf", ylmo%dyn%now%ux_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in x", &
            standard_name="land_ice_surface_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelsurf", ylmo%dyn%now%uy_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in y", &
            standard_name="land_ice_surface_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "zvelsurf", uz_s_masked, &
            units="m s-1", long_name="Surface velocity in z", &
            standard_name="land_ice_surface_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelbase", ylmo%dyn%now%ux_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in x", &
            standard_name="land_ice_basal_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelbase", ylmo%dyn%now%uy_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in y", &
            standard_name="land_ice_basal_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "zvelbase", uz_b_masked, &
            units="m s-1", long_name="Basal velocity in z", &
            standard_name="land_ice_basal_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelmean", ux_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in x", &
            standard_name="land_ice_vertical_mean_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelmean", uy_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in y", &
            standard_name="land_ice_vertical_mean_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litemptop", T_top_ice, &
            units="K", long_name="Surface temperature", &
            standard_name="temperature_at_top_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempavg", T_avg, &
            units="K", long_name="Depth-averaged ice temperature", &
            standard_name="land_ice_temperature", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempbotgr", T_base_grnd, &
            units="K", long_name="Basal temperature beneath grounded ice sheet", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempbotfl", T_base_flt, &
            units="K", long_name="Basal temperature beneath floating ice shelf", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempgradgr", dTdz_base_grnd, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath grounded ice sheet", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempgradfl", dTdz_base_flt, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath floating ice shelf", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "strbasemag", ylmo%dyn%now%taub, &
            units="Pa", long_name="Basal drag", &
            standard_name="land_ice_basal_drag", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftgif", ylmo%tpo%now%f_ice, &
            units="1", long_name="Land ice area fraction", &
            standard_name="land_ice_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftgrf", ylmo%tpo%now%f_grnd, &
            units="1", long_name="Grounded ice sheet area fraction", &
            standard_name="grounded_ice_sheet_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftflf", MAX(ylmo%tpo%now%f_ice - ylmo%tpo%now%f_grnd, 0.0_wp), &
            units="1", long_name="Floating ice sheet area fraction", &
            standard_name="floating_ice_shelf_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ====================================================================
        ! 2-D FL variables  (yearly-average flux)
        ! ====================================================================
        
        ! Time-invariant in most setups; written without time dimension.
        call nc_write(filename, "hfgeoubed", ylmo%bnd%Q_geo * 1.0e3_wp, &
            units="W m-2", long_name="Geothermal heat flux", &
            standard_name="upward_geothermal_heat_flux_in_land_ice", &
            dim1="xc", dim2="yc", start=[1,1], ncid=ncid)
        
        call nc_write(filename, "acabf", ylmo%tpo%now%smb * esm_correction, &
            units="kg m-2 s-1", long_name="Surface mass balance flux", &
            standard_name="land_ice_surface_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "libmassbfgr", bmb_grnd_masked * esm_correction, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath grounded ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "libmassbffl", bmb_shlf_masked * esm_correction, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath floating ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "dlithkdt", ylmo%tpo%now%dHidt / yr_to_sec, &
            units="m s-1", long_name="Ice thickness imbalance", &
            standard_name="tendency_of_land_ice_thickness", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! Kinematic flux through ice-front cells (mask_frnt==1, floating).
        ! Uses same formula as CalvingMIP
        call nc_write(filename, "licalvf", flux_clv_2d, &
            units="kg m-2 s-1", long_name="Calving flux", &
            standard_name="land_ice_specific_mass_flux_due_to_calving", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ligroundf : Grounding-line flux                          [MANDATORY]
        call nc_write(filename, "ligroundf", flux_grl_2d, &
            units="kg m-2 s-1", long_name="Grounding line flux", &
            standard_name="land_ice_specific_grounding_line_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! tfbase : Thermal forcing at ice base, floating           [optional]
        call nc_write(filename, "tfbase", tfbase, &
            units="K", long_name="Thermal forcing at the ice base", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_2D_cmip
        
        
    subroutine write_step_1D_cmip(ylmo, mshlf, filename, time)
        ! Writes all mandatory scalar ISMIP7 (ISM_2026) variables.
        ! Flux computation follows yelmo_calving.f90: kinematic uxy_bar*H_ice*rho_ice.
        ! -------------------------------------------------------------------------
        
        implicit none
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        character(len=*),     intent(IN) :: filename
        real(wp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        type(yregions_class) :: reg
        
        integer  :: ncid, n
        real(wp) :: rho_ice, density_corr, m3yr_to_kgs, esm_correction, yr_to_sec
        real(wp) :: dx, dy
        
        real(wp) :: smb_tot          ! total SMB           [m3 yr-1]
        real(wp) :: bmb_grnd_tot     ! total BMB grounded  [m3 yr-1]
        real(wp) :: bmb_shlf_t       ! total BMB floating  [m3 yr-1]
        real(wp) :: flux_grl         ! total GL flux       [kg yr-1]
        real(wp) :: flux_clv         ! total calving flux  [kg yr-1]
        
        logical, allocatable :: mask_tot(:,:)
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:)
        logical, allocatable :: mask_grl(:,:)    ! grounding-line cells
        logical, allocatable :: mask_frnt(:,:)   ! ice-front cells
        
        ! ---- allocate -------------------------------------------------------
        allocate(mask_tot  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grnd (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_flt  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grl  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_frnt (ylmo%grd%G%nx, ylmo%grd%G%ny))
        
        ! ---- unit conversions -----------------------------------------------
        rho_ice        = 917.0_wp
        m3yr_to_kgs    = 3.2e-5_wp
        density_corr   = rho_ice / 1000.0_wp
        esm_correction = m3yr_to_kgs * density_corr
        yr_to_sec      = 31556952.0_wp
        
        dx = ylmo%grd%G%dx
        dy = ylmo%grd%G%dy
        
        ! ---- masks (updated to match yelmo_calving.f90) ---------------------
        mask_tot  = (ylmo%tpo%now%H_ice .gt. 0.0_wp)
        mask_grnd = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp)
        mask_flt  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
        
        ! Grounding-line cells: grounded ice where mask_grz == 0
        mask_grl  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp &
                        .and. ylmo%tpo%now%mask_grz .eq. 0.0_wp)
        
        ! Ice-front cells: floating ice where mask_frnt == 1
        mask_frnt = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp &
                        .and. ylmo%tpo%now%mask_frnt .eq. 1.0_wp)
        
        ! ---- regional object (pre-computed by Yelmo) ------------------------
        reg = ylmo%reg
        
        ! ---- integrated fluxes ----------------------------------------------
        ! SMB and BMB: area-integrated [m3 yr-1], converted via esm_correction
        smb_tot      = sum(ylmo%bnd%smb,       mask=mask_tot)  * (dx * dy)
        bmb_grnd_tot = sum(ylmo%tpo%now%bmb,   mask=mask_grnd) * (dx * dy)
        bmb_shlf_t   = sum(ylmo%tpo%now%bmb,   mask=mask_flt)  * (dx * dy)
        
        ! Grounding-line flux: kinematic, [kg yr-1]
        ! Matches yelmo_calving.f90: uxy_bar * H_ice * rho_ice * dx
        if (count(mask_grl) .gt. 0) then
            flux_grl = sum(ylmo%dyn%now%uxy_bar * ylmo%tpo%now%H_ice * rho_ice, mask=mask_grl) * dx
        else
            flux_grl = 0.0_wp
        end if
        
        ! Calving flux: sum of cmb_flt and cmb_grnd, [kg yr-1]
        flux_clv = sum(ylmo%tpo%now%cmb_flt + ylmo%tpo%now%cmb_grnd) * (dx * dy)  ! [m3 yr-1]
        
        ! ---- open file & find time index ------------------------------------
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! ====================================================================
        ! Scalar ST variables  (snapshot)
        ! ====================================================================
        
        ! lim : Total ice mass                                     [MANDATORY]
        call nc_write(filename, "lim", reg%V_ice * rho_ice * 1.0e9_wp, &
            units="kg", long_name="Total ice mass", &
            standard_name="land_ice_mass", &
            dim1="time", start=[n], ncid=ncid)
        
        ! limnsw : Mass above floatation                           [MANDATORY]
        call nc_write(filename, "limnsw", reg%V_sl * rho_ice * 1.0e9_wp, &
            units="kg", long_name="Mass above floatation", &
            standard_name="land_ice_mass_not_displacing_sea_water", &
            dim1="time", start=[n], ncid=ncid)
        
        ! iareagr : Grounded ice area                              [MANDATORY]
        call nc_write(filename, "iareagr", reg%A_ice_g * 1.0e6_wp, &
            units="m2", long_name="Grounded ice area", &
            standard_name="grounded_ice_sheet_area", &
            dim1="time", start=[n], ncid=ncid)
        
        ! iareafl : Floating ice area                              [MANDATORY]
        call nc_write(filename, "iareafl", reg%A_ice_f * 1.0e6_wp, &
            units="m2", long_name="Floating ice area", &
            standard_name="floating_ice_shelf_area", &
            dim1="time", start=[n], ncid=ncid)
        
        ! ====================================================================
        ! Scalar FL variables  (yearly-average flux)
        ! ====================================================================
        
        ! tendacabf : Total SMB flux                               [MANDATORY]
        call nc_write(filename, "tendacabf", smb_tot * esm_correction, &
            units="kg s-1", long_name="Total SMB flux", &
            standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbfgr : Total BMB flux, grounded               [MANDATORY]
        call nc_write(filename, "tendlibmassbfgr", bmb_grnd_tot * esm_correction, &
            units="kg s-1", long_name="Total BMB flux beneath grounded ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbffl : Total BMB flux, floating               [MANDATORY]
        call nc_write(filename, "tendlibmassbffl", bmb_shlf_t * esm_correction, &
            units="kg s-1", long_name="Total BMB flux beneath floating ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlicalvf : Total calving flux                         [MANDATORY]
        call nc_write(filename, "tendlicalvf", flux_clv * esm_correction, &
            units="kg s-1", long_name="Total calving flux", &
            standard_name="tendency_of_land_ice_mass_due_to_calving", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlifmassbf : Total ice-front melt flux                [MANDATORY]
        ! ISMIP7-2026 separates this from calving. The kinematic approach does
        ! TODO: replace with sum(ylmo%tpo%now%fmb * ...) when field is confirmed.
        call nc_write(filename, "tendlifmassbf", 0.0_wp, &
            units="kg s-1", long_name="Total ice front melting flux", &
            standard_name="tendency_of_land_ice_mass_due_to_ice_front_melting", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendligroundf : Total grounding-line flux                [MANDATORY]
        ! Kinematic, convert from kg yr-1 to kg s-1.
        call nc_write(filename, "tendligroundf", flux_grl / yr_to_sec, &
            units="kg s-1", long_name="Total grounding line flux", &
            standard_name="tendency_of_grounded_ice_sheet_mass", &
            dim1="time", start=[n], ncid=ncid)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_1D_cmip

end module yelmox_esm_output
