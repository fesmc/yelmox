module yelmox_isos_c_api
  ! C API around FastIsostasy's isostasy + barystatic-sea-level coupling
  ! (libs/yelmox_domain.f90's step_isostasy/couple_isostasy_to_yelmo), exposing
  ! the classic yelmox_esm driver's own isostasy physics to an external host
  ! (Julia's YelmoMirror) -- same rationale/pattern as yelmox_marshelf_c_api.f90.
  ! Single persistent instance, no multi-alias support.
  !
  ! grid_isos == grid_yelmo is assumed (true for the ANT-32KM ISMIP7 setup this
  ! was built for -- yelmox_domain.f90 defaults grid_isos to grid_yelmo when
  ! unset): the remap steps the real driver does between grid_isos and
  ! grid_yelmo are skipped here, isos%out fields are used directly at
  ! (nx, ny) = the yelmo grid. dwdt_corr (grid-resolution relaxation) is an
  ! optional isos_update arg, taken here directly from yelmo bnd%dzbdt_corr.

  use fastisostasy
  use isostasy_defs, only: wp

  implicit none

  type(isos_class), target :: isos1
  type(bsl_class),  target :: bsl1

contains

  subroutine isos_c_init(filename, group, nx, ny, dx, dy, time_rel_init) &
      bind(C, name="isos_init")
    use iso_c_binding
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(in) :: group(*)
    integer(c_int), value         :: nx, ny
    real(c_double), value         :: dx, dy
    real(c_double), value         :: time_rel_init   ! ts%time_rel (constant in "const" spinups)

    ! Shared driver-owned barystatic sea level -- driver calls bsl_init/update
    ! once at program start, before domain_init; folded in here since this
    ! wrapper only ever serves one isostasy instance.
    call bsl_init(bsl1, trim(c_to_f_string(filename)), real(time_rel_init, wp))
    call bsl_update(bsl1, real(time_rel_init, wp))

    call isos_init(isos1, trim(c_to_f_string(filename)), trim(c_to_f_string(group)), &
                    nx, ny, real(dx, wp), real(dy, wp))

  end subroutine

  subroutine isos_c_init_ref(z_bed_ref, H_ice_ref, nx, ny) bind(C, name="isos_init_ref")
    use iso_c_binding
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: z_bed_ref(nx, ny)
    real(c_double), intent(in) :: H_ice_ref(nx, ny)

    call isos_init_ref(isos1, real(z_bed_ref, wp), real(H_ice_ref, wp))

  end subroutine

  subroutine isos_c_init_state(z_bed, H_ice, time, nx, ny) bind(C, name="isos_init_state")
    use iso_c_binding
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: z_bed(nx, ny)
    real(c_double), intent(in) :: H_ice(nx, ny)
    real(c_double), value      :: time

    call isos_init_state(isos1, real(z_bed, wp), real(H_ice, wp), real(time, wp), bsl1)

  end subroutine

  ! Classic restart branch (domain_startup -> bsl_startup + domain_restart_read): restore the shared
  ! sea level from <bundle>/bsl_restart.nc, then initialise isostasy from <bundle>/isos_restart.nc
  ! (state AND reference are read from the bundle, so no isos_init_ref call is needed).
  subroutine isos_c_init_state_restart(fldr, z_bed, H_ice, time, time_rel, nx, ny) &
      bind(C, name="isos_init_state_restart")
    use iso_c_binding
    character(c_char), intent(in) :: fldr(*)
    integer(c_int), value         :: nx, ny
    real(c_double), intent(in)    :: z_bed(nx, ny)
    real(c_double), intent(in)    :: H_ice(nx, ny)
    real(c_double), value         :: time, time_rel

    character(len=1028) :: f

    f = trim(c_to_f_string(fldr))
    call bsl_restart_read(bsl1, trim(f)//"/bsl_restart.nc")
    call bsl_update(bsl1, real(time_rel, wp))

    isos1%par%use_restart = .true.
    isos1%par%restart     = trim(f)//"/isos_restart.nc"
    call isos_init_state(isos1, real(z_bed, wp), real(H_ice, wp), real(time, wp), bsl1)

  end subroutine

  subroutine isos_c_update(H_ice, dwdt_corr, time, time_rel, nx, ny) bind(C, name="isos_update")
    use iso_c_binding
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: H_ice(nx, ny)
    real(c_double), intent(in) :: dwdt_corr(nx, ny)   ! yelmo bnd%dzbdt_corr
    real(c_double), value      :: time       ! ts%time     -> isostasy
    real(c_double), value      :: time_rel   ! ts%time_rel -> shared sea level

    ! Same order as the driver: bsl_update(ts%time_rel) once per step, then
    ! isos_update(ts%time) with the dzbdt correction.
    call bsl_update(bsl1, real(time_rel, wp))
    call isos_update(isos1, real(H_ice, wp), real(time, wp), bsl1, dwdt_corr=real(dwdt_corr, wp))

  end subroutine

  subroutine isos_get_var2D(v2D, nx, ny, name) bind(C, name="isos_get_var2D")
    use iso_c_binding
    integer(c_int), value         :: nx, ny
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(out)   :: v2D(nx, ny)

    character(len=56) :: f_name
    f_name = trim(c_to_f_string(name))

    select case(trim(f_name))
      case("z_bed"); v2D = real(isos1%out%z_bed, c_double)
      case("z_ss");  v2D = real(isos1%out%z_ss,  c_double)
      case("w");     v2D = real(isos1%out%w,     c_double)
      case("we");    v2D = real(isos1%out%we,    c_double)
      case default
        write(*,*) "isos_get_var2D:: Error: variable not recognized: ", trim(f_name)
        stop
    end select

  end subroutine

  function c_to_f_string(c_str) result(f_str)
    use iso_c_binding
    character(c_char), intent(in) :: c_str(*)
    character(len=1028)           :: f_str
    integer                       :: i

    f_str = " "
    do i = 1, 1028
      if (c_str(i) == c_null_char) exit
      f_str(i:i) = c_str(i)
    end do
  end function

end module yelmox_isos_c_api
