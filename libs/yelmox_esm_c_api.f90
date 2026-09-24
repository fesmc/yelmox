module yelmox_esm_c_api
  ! C API around the classic yelmox_esm driver climate + ocean coupling
  ! (yelmox_esm.f90: esm_ctl_load, step_climate_esm, step_marine_shelf_esm),
  ! specialised to grid_clim == grid_smb == grid_mshlf == grid_yelmo == hub grid
  ! (every remap the driver does between them is then a plain copy, so it is
  ! skipped here). Exposes the SAME compiled physics (esm_forcing.f90 +
  ! marine_shelf.f90) to an external host (Julia YelmoMirror).
  !
  ! Only the direct-SMB path (esm.use_smb = True) is wrapped -- the smbpal path
  ! (use_smb = False) stops with a message. Shares the marine-shelf instance
  ! (mshlf1) with yelmox_marshelf_c_api, so both are built into one library
  ! (libs/build_esm_c_api.sh).

  use marine_shelf
  use esm_forcing, only: esm_forcing_class, esm_forcing_init, esm_forcing_update, &
                         esm_variability_update, esm_clim_update
  use nml, only: nml_read
  use yelmox_marshelf_c_api, only: mshlf1, c_to_f_string

  implicit none

  integer, parameter :: wp = kind(1.0)   ! yelmox working precision (single), as in marine_shelf/esm_forcing

  type(esm_forcing_class), target, save :: esm1

  character(len=56)  :: e_run_step
  character(len=256) :: e_domain, e_grid
  logical  :: e_use_esm, e_use_smb, e_use_var, e_use_proj, e_use_hist
  real(wp) :: e_dtt
  real(wp) :: e_time_ref(2), e_time_hist(2), e_time_proj(2), e_time_esm_ref(2)
  character(len=56) :: e_clim_var
  real(wp), allocatable :: smb_out(:,:), tsrf_out(:,:)

contains

  subroutine esm_c_init(path_par, workdir, domain, grid_name, nx, ny) &
      bind(C, name="esm_init")
    use iso_c_binding
    character(c_char), intent(in) :: path_par(*)
    character(c_char), intent(in) :: workdir(*)
    character(c_char), intent(in) :: domain(*)
    character(c_char), intent(in) :: grid_name(*)
    integer(c_int), value         :: nx, ny

    character(len=1028) :: pp, esm_path_par
    character(len=512)  :: par_file
    character(len=56)   :: experiment, esm_name

    pp       = trim(c_to_f_string(path_par))
    e_domain = trim(c_to_f_string(domain))
    e_grid   = trim(c_to_f_string(grid_name))

    ! Mirrors esm_ctl_load (yelmox_esm.f90).
    call nml_read(trim(pp), "ctrl", "run_step", e_run_step)
    call nml_read(trim(pp), "esm", "par_file",   par_file)
    call nml_read(trim(pp), "esm", "experiment", experiment)
    call nml_read(trim(pp), "esm", "esm_name",   esm_name)
    call nml_read(trim(pp), "esm", "use_esm",    e_use_esm)
    call nml_read(trim(pp), "esm", "use_smb",    e_use_smb)
    call nml_read(trim(pp), "esm", "use_var",    e_use_var)
    call nml_read(trim(pp), "esm", "use_proj",   e_use_proj)
    call nml_read(trim(pp), "esm", "use_hist",   e_use_hist)
    call nml_read(trim(pp), "esm", "lapse",        esm1%lapse)
    call nml_read(trim(pp), "esm", "f_p",          esm1%beta_p)
    call nml_read(trim(pp), "esm", "f_ocn",        esm1%f_ocn)
    call nml_read(trim(pp), "esm", "f_polar",      esm1%f_polar)
    call nml_read(trim(pp), "esm", "dT_threshold", esm1%dT_lim)
    call nml_read(trim(pp), "esm", "grid_src",     esm1%grid_src)

    call nml_read(trim(pp), trim(e_run_step), "dtt",          e_dtt)
    call nml_read(trim(pp), trim(e_run_step), "time_ref",     e_time_ref)
    call nml_read(trim(pp), trim(e_run_step), "time_hist",    e_time_hist)
    call nml_read(trim(pp), trim(e_run_step), "time_proj",    e_time_proj)
    call nml_read(trim(pp), trim(e_run_step), "time_esm_ref", e_time_esm_ref)
    call nml_read(trim(pp), trim(e_run_step), "clim_var",     e_clim_var)

    if (.not. e_use_smb) then
      write(*,*) "esm_init:: Error: only esm.use_smb = True (direct SMB) is wrapped."
      stop
    end if

    esm_path_par = trim(c_to_f_string(workdir))//"/"//trim(par_file)
    call esm_forcing_init(esm1, trim(esm_path_par), trim(e_domain), trim(e_grid), &
                          run_type=trim(e_run_step), gcm=trim(esm_name), &
                          experiment=trim(experiment), use_esm=e_use_esm, &
                          use_smb=e_use_smb, use_var=e_use_var, &
                          use_hist=e_use_hist, use_proj=e_use_proj)

    allocate(smb_out(nx, ny), tsrf_out(nx, ny))
    smb_out  = 0.0_wp
    tsrf_out = 0.0_wp

  end subroutine

  ! step_climate_esm, direct-SMB branch, identity grids.
  subroutine esm_c_step_climate(time, z_srf, H_ice, z_bed, f_grnd, z_sl, basins, pd_z_srf, nx, ny) &
      bind(C, name="esm_step_climate")
    use iso_c_binding
    real(c_double), value      :: time
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: z_srf(nx, ny), H_ice(nx, ny), z_bed(nx, ny)
    real(c_double), intent(in) :: f_grnd(nx, ny), z_sl(nx, ny), basins(nx, ny)
    real(c_double), intent(in) :: pd_z_srf(nx, ny)

    real(wp), allocatable :: z_srf_e(:,:), H_ice_e(:,:), z_bed_e(:,:)
    real(wp), allocatable :: f_grnd_e(:,:), z_sl_e(:,:), basins_e(:,:), pd_zsrf_e(:,:)
    real(wp) :: t

    t = real(time, wp)
    z_srf_e = real(z_srf, wp);   H_ice_e = real(H_ice, wp);   z_bed_e = real(z_bed, wp)
    f_grnd_e = real(f_grnd, wp); z_sl_e = real(z_sl, wp);     basins_e = real(basins, wp)
    pd_zsrf_e = real(pd_z_srf, wp)

    ! Step 1: reference climatology (lapse-rate / precip scaling to z_srf).
    call esm_clim_update(esm1, z_srf_e, t, e_time_ref, e_use_smb, trim(e_domain), trim(e_grid))

    if (mshlf1%par%extrap_shlf) then
      write(*,*) "esm_step_climate:: Error: marine_shelf.extrap_shlf = True is not wrapped."
      stop
    end if

    ! Step 2: anomaly fields.
    call esm_forcing_update(esm1, mshlf1, t, e_use_esm, e_time_ref, e_time_hist, e_time_proj, &
                            e_time_esm_ref, trim(e_domain), H_ice_e, basins_e, z_bed_e, f_grnd_e, &
                            z_sl_e, e_use_smb, use_ref_atm=.false., use_ref_ocn=.false.)

    ! Step 3: variability anomaly.
    call esm_variability_update(esm1, mshlf1, t, e_dtt, trim(e_clim_var), e_time_ref, &
                                H_ice_e, basins_e, z_bed_e, f_grnd_e, z_sl_e, e_use_var, &
                                use_ref_atm=.false., use_ref_ocn=.false.)

    ! Direct SMB (esm.use_smb): elevation-corrected annual SMB [mm w.e./yr] + T_srf.
    smb_out  = esm1%smb_ann + sum(esm1%dsmb, dim=3)/12.0_wp - esm1%dsmbdz*(pd_zsrf_e - z_srf_e)
    tsrf_out = sum(esm1%t2m + esm1%dts + esm1%dts_var, dim=3)/12.0_wp
    where(H_ice_e > 0.0_wp .and. tsrf_out > 273.15_wp) tsrf_out = 273.15_wp

  end subroutine

  ! step_marine_shelf_esm, identity grids: ocean forcing at shelf depth + anomalies
  ! into mshlf1%now%T_shlf/S_shlf, then marshelf_update.
  subroutine esm_c_step_marine(H_ice, z_bed, f_grnd, z_sl, regions, basins, dx, nx, ny) &
      bind(C, name="esm_step_marine")
    use iso_c_binding
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: H_ice(nx, ny), z_bed(nx, ny), f_grnd(nx, ny)
    real(c_double), intent(in) :: z_sl(nx, ny), regions(nx, ny), basins(nx, ny)
    real(c_double), value      :: dx

    real(wp), allocatable :: H(:,:), zb(:,:), fg(:,:), zsl(:,:), reg(:,:), bas(:,:)
    real(wp), allocatable :: T_e(:,:), S_e(:,:)

    H = real(H_ice, wp); zb = real(z_bed, wp); fg = real(f_grnd, wp)
    zsl = real(z_sl, wp); reg = real(regions, wp); bas = real(basins, wp)

    allocate(T_e(nx, ny), S_e(nx, ny))
    call marshelf_interp_shelf(T_e, mshlf1, esm1%to_ref%var(:,:,:,1), H, zb, fg, zsl, -esm1%to_ref%z)
    call marshelf_interp_shelf(S_e, mshlf1, esm1%so_ref%var(:,:,:,1), H, zb, fg, zsl, -esm1%so_ref%z)
    T_e = T_e + esm1%dto + esm1%dto_var
    S_e = S_e + esm1%dso + esm1%dso_var

    mshlf1%now%T_shlf = T_e
    mshlf1%now%S_shlf = S_e

    if (trim(e_domain) == "Greenland") then
      mshlf1%now%dT_shlf = T_e + esm1%dto
      mshlf1%par%tf_method = 2
    end if

    call marshelf_update(mshlf1, H, zb, fg, reg, bas, zsl, real(dx, wp))

  end subroutine

  subroutine esm_c_get_var2D(v2D, nx, ny, name) bind(C, name="esm_get_var2D")
    use iso_c_binding
    integer(c_int), value         :: nx, ny
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(out)   :: v2D(nx, ny)

    character(len=56) :: f_name
    f_name = trim(c_to_f_string(name))

    select case(trim(f_name))
      case("smb");     v2D = real(smb_out,  c_double)   ! [mm w.e./yr], as dom%smb%ann%smb
      case("tsrf");    v2D = real(tsrf_out, c_double)   ! [K], as dom%smb%ann%tsrf
      case("dto");     v2D = real(esm1%dto,     c_double)
      case("dso");     v2D = real(esm1%dso,     c_double)
      case("dto_var"); v2D = real(esm1%dto_var, c_double)
      case("Qd_ann");  v2D = real(esm1%Qd_ann,  c_double)
      case default
        write(*,*) "esm_get_var2D:: Error: variable not recognized: ", trim(f_name)
        stop
    end select

  end subroutine

end module yelmox_esm_c_api
