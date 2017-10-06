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
    MAXP     =  200,             &
    MAXL     = 1000,             &
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
    SOLID    =    8,             & ! material state: solid
    CMN_FILE =    7                ! T-FlowS command file (T-FlowS.cmn)

  real, parameter ::  & 
    HUGE=1.e+30, TINY=1.e-64

!----- Unknown type
  TYPE Unknown          
    real,pointer :: n(:)                ! new value
    real,pointer :: o(:), oo(:)         ! old and older then old
    real,pointer :: C(:), Co(:), Coo(:) ! convective fluxes
    real,pointer :: Do(:), Doo(:)       ! difussive fluxes
    real,pointer :: X(:), Xo(:), Xoo(:) ! surfce sources  
    real,pointer :: mean(:)             ! long time average
    real,pointer :: filt(:)             ! long time average
    real,pointer :: q(:)                ! flux of a variable
    real,pointer :: fluc(:) 
    real         :: URF                 ! under relaxation factor
    real         :: Stol                ! solver tolerance
    real         :: bound(128)          ! boundary values
    real         :: init(128)           ! initial values
    real         :: pro(11024)           ! inlfow profile
    real         :: Sigma               ! sigma
  end TYPE Unknown

end MODULE 
