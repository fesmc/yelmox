module yelmox_marshelf_c_api
  ! C API around libs/marine_shelf.f90 -- exposes the classic yelmox_esm
  ! driver's own marine-melt scheme (bmb_method="quad-nl" etc, including
  ! tf_corr) to an external host (Julia's YelmoMirror), so the ISMIP7 mirror
  ! validation script can use the SAME compiled physics the classic driver
  ! uses instead of having none. Single persistent instance (no multi-alias
  ! support, unlike yelmo_c_api.f90's ylmo1/ylmo2) -- add aliasing later if
  ! more than one marine-shelf state per process is ever needed.
  !
  ! Inputs to marshelf_update (H_ice, z_bed, f_grnd, regions, basins, z_sl)
  ! are all already reachable from Julia via the existing yelmo C API
  ! (tpo_H_ice, bnd_z_bed, dyn_f_grnd, bnd_regions, bnd_basins, bnd_z_sl) --
  ! this module doesn't duplicate getters for them, only for marshelf's own
  ! output/state fields.

  use marine_shelf

  implicit none

  ! marine_shelf.f90 declares wp = sp internally but does not export it; match it here.
  integer, parameter :: wp = kind(1.0)

  type(marshelf_class), target :: mshlf1

contains

  subroutine marshelf_c_init(filename, group, nx, ny, domain, grid_name, &
                              regions, basins, xc, yc, dx) &
      bind(C, name="marshelf_init")
    use iso_c_binding
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(in) :: group(*)
    integer(c_int), value         :: nx, ny
    character(c_char), intent(in) :: domain(*)
    character(c_char), intent(in) :: grid_name(*)
    real(c_double), intent(in)    :: regions(nx, ny)
    real(c_double), intent(in)    :: basins(nx, ny)
    real(c_double), intent(in)    :: xc(nx)
    real(c_double), intent(in)    :: yc(ny)
    real(c_double), value         :: dx

    call marshelf_init(mshlf1, c_to_f_string(filename), c_to_f_string(group), nx, ny, &
                        c_to_f_string(domain), c_to_f_string(grid_name), &
                        real(regions, wp), real(basins, wp), &
                        real(xc, wp), real(yc, wp), real(dx, wp))

  end subroutine

  subroutine marshelf_c_update(H_ice, z_bed, f_grnd, regions, basins, z_sl, dx, nx, ny) &
      bind(C, name="marshelf_update")
    use iso_c_binding
    integer(c_int), value      :: nx, ny
    real(c_double), intent(in) :: H_ice(nx, ny)
    real(c_double), intent(in) :: z_bed(nx, ny)
    real(c_double), intent(in) :: f_grnd(nx, ny)
    real(c_double), intent(in) :: regions(nx, ny)
    real(c_double), intent(in) :: basins(nx, ny)
    real(c_double), intent(in) :: z_sl(nx, ny)
    real(c_double), value      :: dx

    call marshelf_update(mshlf1, real(H_ice, wp), real(z_bed, wp), real(f_grnd, wp), &
                          real(regions, wp), real(basins, wp), real(z_sl, wp), real(dx, wp))

  end subroutine

  subroutine marshelf_get_var2D(v2D, nx, ny, name) bind(C, name="marshelf_get_var2D")
    use iso_c_binding
    integer(c_int), value         :: nx, ny
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(out)   :: v2D(nx, ny)

    character(len=56) :: f_name
    f_name = trim(c_to_f_string(name))

    select case(trim(f_name))
      case("bmb_shlf");      v2D = real(mshlf1%now%bmb_shlf,      c_double)
      case("bmb_ref");       v2D = real(mshlf1%now%bmb_ref,       c_double)
      case("bmb_corr");      v2D = real(mshlf1%now%bmb_corr,      c_double)
      case("T_shlf");        v2D = real(mshlf1%now%T_shlf,        c_double)
      case("dT_shlf");       v2D = real(mshlf1%now%dT_shlf,       c_double)
      case("S_shlf");        v2D = real(mshlf1%now%S_shlf,        c_double)
      case("T_fp_shlf");     v2D = real(mshlf1%now%T_fp_shlf,     c_double)
      case("tf_shlf");       v2D = real(mshlf1%now%tf_shlf,       c_double)
      case("tf_corr");       v2D = real(mshlf1%now%tf_corr,       c_double)
      case("tf_corr_basin"); v2D = real(mshlf1%now%tf_corr_basin, c_double)
      case("mask_ocn");      v2D = real(mshlf1%now%mask_ocn,      c_double)
      case default
        write(*,*) "marshelf_get_var2D:: Error: variable not recognized: ", trim(f_name)
        stop
    end select

  end subroutine

  subroutine marshelf_set_var2D(v2D, nx, ny, name) bind(C, name="marshelf_set_var2D")
    use iso_c_binding
    integer(c_int), value         :: nx, ny
    character(c_char), intent(in) :: name(*)
    real(c_double), intent(in)    :: v2D(nx, ny)

    character(len=56) :: f_name
    f_name = trim(c_to_f_string(name))

    select case(trim(f_name))
      ! tf_corr is the intended write target -- an external host (the
      ! optimize_tf_corr! port) drives it, marshelf_update then folds it
      ! into bmb_shlf via whatever bmb_method is configured (matches how
      ! the classic driver's own tf_corr feeds this same field).
      case("tf_corr"); mshlf1%now%tf_corr = real(v2D, wp)
      case default
        write(*,*) "marshelf_set_var2D:: Error: variable not recognized or not settable: ", trim(f_name)
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

end module yelmox_marshelf_c_api
