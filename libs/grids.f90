
module grids
     
    use coord 
    
    implicit none

contains

    subroutine grid_select(grid,grid_name)

        implicit none 

        type(grid_class), intent(OUT) :: grid  
        character(len=*), intent(IN)  :: grid_name 

        select case(trim(grid_name))

            ! GREENLAND DOMAINS =======================

            case("ESPG-3413-20KM")
                call grid_init(grid,name="ESPG-3413-20KM",mtype="polar_stereographic", &
                              units="kilometers",lon180=.TRUE., &
                              x0=-720.d0,dx=20.0d0,nx=85,y0=-3450.d0,dy=20.0d0,ny=145, &
                              lambda=-45.d0,phi=70.d0,alpha=20.0d0)
            
            case("ESPG-3413-10KM")
                call grid_init(grid,name="ESPG-3413-10KM",mtype="polar_stereographic", &
                              units="kilometers",lon180=.TRUE., &
                              x0=-720.d0,dx=10.0d0,nx=169,y0=-3450.d0,dy=10.0d0,ny=289, &
                              lambda=-45.d0,phi=70.d0,alpha=20.0d0)
            
            case("ESPG-3413-5KM")
                call grid_init(grid,name="ESPG-3413-5KM",mtype="polar_stereographic", &
                              units="kilometers",lon180=.TRUE., &
                              x0=-720.d0,dx=5.0d0,nx=337,y0=-3450.d0,dy=5.0d0,ny=577, &
                              lambda=-45.d0,phi=70.d0,alpha=20.0d0)
            
            case("ESPG-3413-2KM")
                call grid_init(grid,name="ESPG-3413-2KM",mtype="polar_stereographic", &
                              units="kilometers",lon180=.TRUE., &
                              x0=-720.d0,dx=2.0d0,nx=841,y0=-3450.d0,dy=2.0d0,ny=1441, &
                              lambda=-45.d0,phi=70.d0,alpha=20.0d0)
            
            case("ESPG-3413-1KM")
                call grid_init(grid,name="ESPG-3413-1KM",mtype="polar_stereographic", &
                              units="kilometers",lon180=.TRUE., &
                              x0=-720.d0,dx=1.0d0,nx=1681,y0=-3450.d0,dy=1.0d0,ny=2881, &
                              lambda=-45.d0,phi=70.d0,alpha=20.0d0)
            
            case("GRL-120KM")
                call grid_init(grid,name="GRL-120KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=120.d0,nx=16,dy=120.d0,ny=26, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-40KM")
                call grid_init(grid,name="GRL-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=46,dy=40.d0,ny=76, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-20KM")
                call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-10KM")
                call grid_init(grid,name="GRL-10KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=10.d0,nx=181,dy=10.d0,ny=301, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-5KM")
                call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)
            
            case("GRL-2KM")
                call grid_init(grid,name="GRL-2KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=2.d0,nx=901,dy=2.d0,ny=1501, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("GRL-1KM")
                call grid_init(grid,name="GRL-1KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=1.d0,nx=1801,dy=1.d0,ny=3001, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

            case("Bamber01-20KM")
                call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
                               lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                               lambda=-39.d0,phi=90.d0,alpha=19.d0)
    
            case("Bamber01-10KM")
                call grid_init(grid,name="Bamber01-10KM",mtype="polar_stereographic",units="kilometers", &
                               lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
                               lambda=-39.d0,phi=90.d0,alpha=19.d0)

            ! ANTARCTICA DOMAINS ======================= 

            case("ANT-80KM")
                call grid_init(grid,name="ANT-80KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=80.d0,nx=79,dy=80.d0,ny=74,lambda=0.d0,phi=-71.d0)

            case("ANT-40KM")
                call grid_init(grid,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

            case("ANT-20KM")
                call grid_init(grid,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=20.d0,nx=313,dy=20.d0,ny=293,lambda=0.d0,phi=-71.d0)

            case("ANT-10KM")
                call grid_init(grid,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

            case("ANT-5KM")
                call grid_init(grid,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=5.d0,nx=1249,dy=5.d0,ny=1169,lambda=0.d0,phi=-71.d0)

            case("ANT-1KM")
                call grid_init(grid,name="ANT-1KM",mtype="polar_stereographic",units="kilometers", &
                       lon180=.TRUE.,dx=1.d0,nx=6241,dy=1.d0,ny=5841,lambda=0.d0,phi=-71.d0)

            ! NORTH DOMAINS ======================= 

            case("NH-40KM")
                call grid_init(grid,name="NH-40KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=40.d0,nx=225,dy=40.d0,ny=211, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-20KM")
                call grid_init(grid,name="NH-20KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=20.d0,nx=449,dy=20.d0,ny=421, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-10KM")
                call grid_init(grid,name="NH-10KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=10.d0,nx=897,dy=10.d0,ny=841, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-5KM")
                call grid_init(grid,name="NH-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=1793,dy=5.d0,ny=1681, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case("NH-1KM")
                call grid_init(grid,name="NH-1KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=1.d0,nx=8961,dy=1.d0,ny=8401, &
                               lambda=-53.d0,phi=78.d0,alpha=32.7d0)

            case DEFAULT
                write(*,*) "grid_select:: error: grid name not recognized: "//trim(grid_name)
                stop 

        end select

        return 

    end subroutine grid_select
    
end module grids 
