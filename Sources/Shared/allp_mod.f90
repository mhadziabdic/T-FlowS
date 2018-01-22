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

  !---------------------------------!
  !   Twp handy logical constants   !
  !---------------------------------!
  integer, parameter :: YES =  0
  integer, parameter :: NO  = -1

  integer, parameter ::          & 
    FLUID    =    7,             & ! material state: fluid
    SOLID    =    8,             & ! material state: solid
    CMN_FILE =    7                ! T-FlowS command file (T-FlowS.cmn)

  !----------------------------------------!
  !   A few handy mathematical constants   !
  !----------------------------------------!
  real, parameter :: HUGE       = 1.e+30
  real, parameter :: TINY       = 1.e-30
  real, parameter :: PI         = 3.14159265359
  real, parameter :: ONE_THIRD  = 0.33333333333333333
  real, parameter :: TWO_THIRDS = 1.0 - ONE_THIRD

end module 
