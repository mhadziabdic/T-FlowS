!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!      Parameter definitions      !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE allp_mod

  implicit none

  integer, parameter ::          & 
    MAXP =   200, MAXL = 1000,   &
    MAXPRO   = 1024,             & ! max. n. of processors    
    INFLOW   =    1,             & ! boundary condition       
    WALL     =    2,             & ! boundary condition
    OUTFLOW  =    3,             & ! boundary condition
    SYMMETRY =    4,             & ! boundary condition
    CONVECT  =    5,             & ! boundary condition
    PRESSURE = 2312,             & ! boundary condition 
    PERIODIC = 2313,             & ! boundary condition
    BUFFER   = 2311,             & ! boundary condition 
                                   ! BUFFER CANNOT BE 11 AS THIS 
                                   ! MAKES PROBLEM IF CASE HAS MORE THAN 10 BOUNDARY
                                   ! CONDITIONS. 11th BOUNDARY CONDITION IS CONFUSED
                                   ! FOR BUFFER IN CASE BUFFER = 11. 
    WALLFL   =    6                ! boundary condition

  integer, parameter ::          & 
    FLUID    =    7,             & ! material state: fluid
    SOLID    =    8                ! material state: solid

  real, parameter ::  & 
    HUGE=1.e+30, TINY=1.e-64

!----- Unknown type
  TYPE Unknown          
    real,POINTER :: n(:)                ! new value
    real,POINTER :: o(:), oo(:)         ! old and older then old
    real,POINTER :: C(:), Co(:), Coo(:) ! convective fluxes
    real,POINTER :: Do(:), Doo(:)       ! difussive fluxes
    real,POINTER :: X(:), Xo(:), Xoo(:) ! surfce sources  
    real,POINTER :: mean(:)             ! long time average
    real,POINTER :: filt(:)             ! long time average
    real,POINTER :: q(:)                ! flux of a variable
    real,POINTER :: fluc(:) 
    real         :: URF                 ! under relaxation factor
    real         :: Stol                ! solver tolerance
    real         :: bound(128)          ! boundary values
    real         :: init(128)           ! initial values
    real         :: pro(11024)           ! inlfow profile
    real         :: Sigma               ! sigma
  end TYPE Unknown

end MODULE 
