module pico_geometry
    ! Module to simulate the geometry of ice shelves via PICO

    use nml 
    use ncio 

!     use yelmo_defs, only : sp, dp, prec, rho_ice, rho_w, rho_sw, g, parse_path 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    private

    public :: pico_calc_shelf_extent
    public :: pico_calc_area_boxes
    public :: pico_calc_shelf_boxes
    public :: calc_grline
    public :: calc_margin

contains 

    ! Compute shelf extent
    subroutine pico_calc_shelf_extent(d_shlf,d_if,r_shlf,H_ice,f_grnd,is_grline,is_margin,dx)

        implicit none

        real(prec), intent(INOUT) :: d_shlf(:,:),d_if(:,:),r_shlf(:,:)
        real(prec), intent(IN) :: H_ice(:,:),f_grnd(:,:)
        logical :: is_grline(:,:),is_margin(:,:)
        real(prec), intent(IN) :: dx

        ! Set every step to zero to initialize
        d_shlf = eikonal_equation(is_grline,f_grnd,H_ice)*dx*0.001 ! to km
        d_if = eikonal_equation(is_margin,f_grnd,H_ice)*dx*0.001   ! to km
        r_shlf = 0.0
        where (d_shlf .ne. 0.0 .or. d_if .ne. 0.0) r_shlf =d_shlf/(d_shlf+d_if)

        return

    end subroutine pico_calc_shelf_extent

    ! Function to compute the area of every box
    function pico_calc_area_boxes(boxes,f_grnd,n_box,dx) result(A_box)

        implicit none

        real(prec), intent(IN) :: boxes(:,:),f_grnd(:,:)
        real(prec), intent(IN)  :: n_box, dx
        real(prec) :: A_box(size(boxes,1),size(boxes,2))

       ! Local variables
       real(prec) :: n_count
       integer :: n,i,j

       ! Initialize variable
       A_box = 0.0

       do n=1,int(n_box)
           n_count = 0.0
           do i=1, size(boxes,1)
               do j=1, size(boxes,2)
                   if (boxes(i,j) .eq. n) n_count = n_count + 1.0 -f_grnd(i,j)
               end do
           end do
           ! test count
           where(boxes .eq. n) A_box = n_count*dx*dx ! Esta en m
       end do

       return

    end function pico_calc_area_boxes

    ! Function to compute boxes
    function pico_calc_shelf_boxes(H_ice,f_grnd,d_shlf,r_shlf,n_box,d_max,d_max_basin) result(boxes)

        implicit none

        real(prec), intent(IN) :: H_ice(:,:),f_grnd(:,:),d_shlf(:,:),r_shlf(:,:)
        real(prec), intent(IN) :: n_box,d_max,d_max_basin
        real(prec) :: boxes(size(f_grnd,1),size(f_grnd,2))

        ! Local variables
        integer :: i, j, k
        real(prec) :: n_shlf

        ! Initialize
        n_shlf = 0.0
        boxes  = 0.0

        ! Assign number of boxes
        n_shlf = 1.0 + NINT(SQRT(d_max_basin/d_max)*(n_box-1))

        ! divide ice shelves
        do k=1,int(n_shlf)
            do i=1, size(f_grnd,1)
            do j=1, size(f_grnd,2)
                ! Floating points
                if(f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) then
                    if((r_shlf(i,j) .ge. (1-SQRT((n_shlf-k+1)/n_shlf))) .and. (r_shlf(i,j) .le. (1-SQRT((n_shlf-k)/n_shlf)))) then
                        boxes(i,j) = k
                    end if
                end if
            end do
            end do
        end do

        return

    end function pico_calc_shelf_boxes   
 
    ! ======================================
    ! 
    ! Distance functions + GRL definition 
    ! 
    ! ======================================
 
    function calc_grline(f_grnd,H_ice) result(is_grline)
        ! Determine the grounding line given the grounded fraction f_grnd
        ! ie, is_grline is true for a floating point or partially floating  
        ! point with grounded neighbors

        implicit none 

        real(prec), intent(IN) :: f_grnd(:,:),H_ice(:,:)
        logical :: is_grline(size(f_grnd,1),size(f_grnd,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(is_grline,1)
        ny = size(is_grline,2)

        is_grline = .FALSE. 
        do i = 2, nx-1 
        do j = 2, ny-1 
          
            ! jablasco: grl definition in PICO first purely floating point ->
            ! avoid melting in grounding line (box1) 
            ! Floating point or partially floating point with grounded neighbors
            if ((H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .eq. 0.0) .and. &
                (f_grnd(i-1,j) .gt. 0.0 .or. f_grnd(i+1,j) .gt. 0.0 .or. &
                 f_grnd(i,j-1) .gt. 0.0 .or. f_grnd(i,j+1) .gt. 0.0) ) then 
                is_grline(i,j) = .TRUE. 

            end if 
            
        end do 
        end do 

        return 

    end function calc_grline

    function calc_margin(H_ice) result(is_margin)
        ! Determine the ice margin

        implicit none

        real(prec), intent(IN) :: H_ice(:,:)
        logical :: is_margin(size(H_ice,1),size(H_ice,2))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(is_margin,1)
        ny = size(is_margin,2)

        is_margin = .FALSE.
        do i = 2, nx-1
        do j = 2, ny-1

            ! Floating point or partially floating point with ocean neighbors
            if (H_ice(i,j) .gt. 0.0 .and. &
                (H_ice(i-1,j) .eq. 0.0 .or. H_ice(i+1,j) .eq. 0.0 .or. &
                 H_ice(i,j-1) .eq. 0.0 .or. H_ice(i,j+1) .eq. 0.0) ) then
                is_margin(i,j) = .TRUE.

            end if

        end do
        end do

        return

    end function calc_margin

    function eikonal_equation(mask,f_grnd,H_ice) result(dists)

        implicit none

        ! Logical mask, either is_grline or is_margin
        logical,    intent(IN) :: mask(:,:)
        real(prec), intent(IN) :: f_grnd(:,:), H_ice(:,:)
        real(prec) :: dists(size(mask,1),size(mask,2))

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: loop, current_label
 
        nx = size(mask,1)
        ny = size(mask,2)

        loop = 1.0
        current_label = 1.0
        ! Assign to dists mask
        dists = 0.0
        where(mask) dists = 1.0

        do while(loop .ne. 0.0)
            loop = 0.0

            do i = 2, nx-1
            do j = 2, ny-1

                ! Only floating points with no distance (no grounding line)
                if ((f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) .and. dists(i,j) .eq. 0.0) then
                    if(((dists(i-1,j) .eq. current_label)) .or. (dists(i+1,j) .eq. current_label) .or. (dists(i,j-1) .eq. current_label) .or. (dists(i,j+1) .eq. current_label)) then
                        dists(i,j) = current_label + 1.0
                        loop = 1.0
                    end if
                end if
            end do
            end do

            current_label = current_label+1.0

        end do

        ! jablasco: correct distance! (distance grl is 0 not 1)
        where(dists .gt. 0.0) dists = dists - 1.0

        return

    end function eikonal_equation

end module pico_geometry
