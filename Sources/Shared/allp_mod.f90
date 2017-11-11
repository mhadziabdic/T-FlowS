!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!      Parameter definitions      !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module allp_mod

  implicit none

  integer, parameter :: YES =  0
  integer, parameter :: NO  = -1

  integer, parameter ::          & 
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

  real, parameter ::                   & 
    HUGE       = 1.e+30,               &
    TINY       = 1.e-30,               &
    PI         = 3.14159265359,        &
    ONE_THIRD  = 0.33333333333333333,  &
    TWO_THIRDS = 1.0 - ONE_THIRD

end module 
