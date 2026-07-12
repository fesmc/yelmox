program test_snapclim_ref
    ! Validation reference dumper for the snapclim -> snapclim2 refactor.
    !
    ! Drives the legacy `snapclim` on a real config (Greenland GRL-16KM,
    ! atm_type=snap_1ind_new) with synthetic-but-deterministic z_srf/basins, and
    ! writes now%{...} to NetCDF at a few times. snapclim reads the REAL snapshot
    ! data from ice_data; only z_srf/basins are supplied here so that the identical
    ! inputs can later be fed to snapclim2 for a field-by-field diff.
    !
    ! Run from the yelmox root so relative paths (yelmox/, ice_data/, input/, logs/)
    ! resolve.

    use snapclim

    implicit none

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: nx = 106, ny = 181

    type(snapclim_class) :: snp
    real(sp) :: z_srf(nx,ny), basins(nx,ny)
    real(sp) :: xc(nx), yc(ny)
    real(sp) :: times(3)
    character(len=256) :: path_par, outfile, domain, grid
    integer :: i, j, it, nt

    path_par = "yelmox/yelmox_Greenland.nml"
    domain   = "Greenland"
    grid     = "GRL-16KM"
    outfile  = "logs/snapclim_ref.nc"

    ! Synthetic, deterministic surface elevation (a smooth dome) to exercise the
    ! lapse-rate correction, and uniform basins.
    do j = 1, ny
    do i = 1, nx
        z_srf(i,j) = 2500.0*exp(-( ((i-nx/2.0)/40.0)**2 + ((j-ny/2.0)/60.0)**2 ))
    end do
    end do
    basins = 1.0

    do i = 1, nx; xc(i) = (i-1)*16.0; end do
    do j = 1, ny; yc(j) = (j-1)*16.0; end do

    times = [0.0, -21000.0, -120000.0]
    nt = size(times)

    call snapclim_init(snp, trim(path_par), trim(domain), trim(grid), nx, ny, basins, group="snap")

    call snapclim_write_init(snp, trim(outfile), xc, yc, times(1), units="years")

    do it = 1, nt
        call snapclim_update(snp, z_srf=z_srf, time=times(it), domain=trim(domain), basins=basins)
        call snapclim_write_step(snp, trim(outfile), times(it))
        write(*,"(a,f12.1,a,f10.4,a,f10.4)") " time=", times(it), &
            "  ta_ann(mid)=", snp%now%ta_ann(nx/2,ny/2), &
            "  pr_ann(mid)=", snp%now%pr_ann(nx/2,ny/2)
    end do

    write(*,*) "wrote ", trim(outfile)

end program test_snapclim_ref
