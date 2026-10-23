
module interp_time
    
    use interp1D 
    
    implicit none 

    private
    public :: convert_monthly_daily_1D, convert_monthly_daily_3D

contains

    subroutine convert_monthly_daily_3D(mon3D,day3D,days)

        implicit none 

        double precision, dimension(:,:,:) :: mon3D, day3D 
        double precision, dimension(:,:), allocatable :: mon0, mon13
        integer, dimension(:), optional :: days 
        integer :: nx, ny, nm, nd, nk 
        integer :: i, j, k 

        nx = size(mon3D,1)
        ny = size(mon3D,2)
        nm = size(mon3D,3)

        nd = size(day3D,3)

        if (allocated(mon0))  deallocate(mon0)
        if (allocated(mon13)) deallocate(mon13) 
        allocate(mon0(nx,ny),mon13(nx,ny))

        mon0  = mon3D(:,:,12)
        mon13 = mon3D(:,:,1)

        !$omp parallel do collapse(2) private(i,j)
        do j = 1, ny 
        do i = 1, nx 
            day3D(i,j,:) = convert_monthly_daily_1D(mon3D(i,j,:), &
                             mon0(i,j),mon13(i,j),nd,days)
        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine convert_monthly_daily_3D

    function convert_monthly_daily_1D(dat_mon,dat_m0,dat_m13,nd,days) result(dat_day)
        
        implicit none

        double precision, dimension(:), intent(IN)   :: dat_mon
        double precision,               intent(IN)   :: dat_m0, dat_m13
        integer :: nd 
        integer,          dimension(:), intent(IN), optional :: days 
        double precision, dimension(nd)  :: dat_day
        double precision, dimension(:), allocatable  :: dat_mon14
        double precision, dimension(:), allocatable :: x, xout 
        integer :: i 

        allocate(x(14),xout(nd)) 
        allocate(dat_mon14(14))

        ! Define x-values for converting to daily data from months using splines
        do i = 1, 14
            x(i) = dble(i-1)*30-15
        end do 
        
        if (present(days)) then 
            xout = days 
        else 
            do i = 1, nd
                xout(i) = dble(i)
            end do 
        end if 

        if (.not. present(days) .and. nd .ne. 360) then 
            write(*,*) "convert_monthly_daily_1D:: ", &
            "Error: This routine only works with 360-day years &
             &unless specific days are provided."
            stop 
        else if (nd .ne. size(xout)) then 
            write(*,*) "convert_monthly_daily_1D:: ", &
            "Error: Number of days given and size of fields do not match."
            stop 
        end if 

        if (minval(xout) .lt. x(1) .or. maxval(xout) .gt. x(14)) then 
            write(*,*) "convert_monthly_days_1D:: ", &
            "Error: currently this routine only works with 360-day years."
            stop 
        end if 

        dat_mon14(1)    = dat_m0
        dat_mon14(14)   = dat_m13 
        dat_mon14(2:13) = dat_mon 

        
        dat_day = interp_spline(x,dat_mon14,xout)

        return 

    end function convert_monthly_daily_1D

end module interp_time 
