program test_snapesm_ref
    ! Validation harness for the snapclim -> snapesm port. Mirrors
    ! tests/test_snapclim_ref.f90 exactly (same GRL-16KM config, same synthetic
    ! z_srf dome / uniform basins, same times), driving `snapesm` instead of
    ! snapclim and dumping now%{...} to logs/snapesm_ref.nc for a field-by-field
    ! diff (tests/diff_nc.py) against logs/snapclim_ref.nc.
    !
    ! Run from the yelmox root so relative paths (ice_data/, input/, logs/, and
    ! the .claude/worktrees/snapclim2/input config) resolve.

    use precision, only : wp, sp
    use ncio
    use snapesm

    implicit none

    integer, parameter :: nx = 106, ny = 181

    type(snapesm_class) :: sc
    real(sp) :: z_srf_sp(nx,ny), basins_sp(nx,ny)
    real(wp) :: z_srf(nx,ny), basins(nx,ny)
    real(wp) :: xc(nx), yc(ny)
    real(wp) :: times(3)
    character(len=256) :: path_par, outfile, domain, grid
    integer :: i, j, it, nt

    path_par = ".claude/worktrees/snapclim2/input/greenland_snp.nml"
    domain   = "Greenland"
    grid     = "GRL-16KM"
    outfile  = "logs/snapesm_ref.nc"

    ! Identical synthetic z_srf dome & uniform basins as test_snapclim_ref (computed
    ! in single precision, then promoted, so the inputs are bit-identical).
    do j = 1, ny
    do i = 1, nx
        z_srf_sp(i,j) = 2500.0*exp(-( ((i-nx/2.0)/40.0)**2 + ((j-ny/2.0)/60.0)**2 ))
    end do
    end do
    basins_sp = 1.0
    z_srf  = real(z_srf_sp, wp)
    basins = real(basins_sp, wp)

    do i = 1, nx; xc(i) = (i-1)*16.0_wp; end do
    do j = 1, ny; yc(j) = (j-1)*16.0_wp; end do

    times = [0.0_wp, -21000.0_wp, -120000.0_wp]
    nt = size(times)

    call snapesm_init(sc, trim(path_par), trim(domain), trim(grid), nx, ny, &
                      times(1), basins, group="snp")

    call dump_init(trim(outfile), sc, xc, yc, times(1))

    do it = 1, nt
        call snapesm_update(sc, z_srf=z_srf, time=times(it), domain=trim(domain), basins=basins)
        call dump_step(trim(outfile), sc, times(it), it)
        write(*,"(a,f12.1,a,f10.4,a,f10.4)") " time=", times(it), &
            "  ta_ann(mid)=", sc%now%ta_ann(nx/2,ny/2), &
            "  pr_ann(mid)=", sc%now%pr_ann(nx/2,ny/2)
    end do

    write(*,*) "wrote ", trim(outfile)

contains

    subroutine dump_init(filename, sc, xc, yc, time0)
        character(len=*),      intent(IN) :: filename
        type(snapesm_class),   intent(IN) :: sc
        real(wp),              intent(IN) :: xc(:), yc(:), time0
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=xc, units="km")
        call nc_write_dim(filename,"yc",    x=yc, units="km")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12, units="month")
        call nc_write_dim(filename,"depth", x=sc%now%depth, units="meters")
        call nc_write_dim(filename,"time",  x=time0,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)
        return
    end subroutine dump_init

    subroutine dump_step(filename, sc, time, n)
        character(len=*),      intent(IN) :: filename
        type(snapesm_class),   intent(IN) :: sc
        real(wp),              intent(IN) :: time
        integer,               intent(IN) :: n
        integer :: ncid
        call nc_open(filename, ncid, writable=.TRUE.)
        call nc_write(filename,"time", time, dim1="time", start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"mask",  sc%now%mask,  dim1="xc",dim2="yc",dim3="time", &
                        start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"z_srf", sc%now%z_srf, dim1="xc",dim2="yc",dim3="time", &
                        start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"tsl", sc%now%tsl, dim1="xc",dim2="yc",dim3="month",dim4="time",&
                        start=[1,1,1,n],count=[nx,ny,12,1],ncid=ncid)
        call nc_write(filename,"tas", sc%now%tas, dim1="xc",dim2="yc",dim3="month",dim4="time",&
                        start=[1,1,1,n],count=[nx,ny,12,1],ncid=ncid)
        call nc_write(filename,"pr",  sc%now%pr,  dim1="xc",dim2="yc",dim3="month",dim4="time",&
                        start=[1,1,1,n],count=[nx,ny,12,1],ncid=ncid)
        call nc_write(filename,"ta_ann", sc%now%ta_ann, dim1="xc",dim2="yc",dim3="time", &
                        start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"ta_sum", sc%now%ta_sum, dim1="xc",dim2="yc",dim3="time", &
                        start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"pr_ann", sc%now%pr_ann, dim1="xc",dim2="yc",dim3="time", &
                        start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"to_ann", sc%now%to_ann, dim1="xc",dim2="yc",dim3="depth",dim4="time",&
                        start=[1,1,1,n],count=[nx,ny,23,1],ncid=ncid)
        call nc_write(filename,"so_ann", sc%now%so_ann, dim1="xc",dim2="yc",dim3="depth",dim4="time",&
                        start=[1,1,1,n],count=[nx,ny,23,1],ncid=ncid)
        call nc_close(ncid)
        return
    end subroutine dump_step

end program test_snapesm_ref
