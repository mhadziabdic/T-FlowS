!======================================================================!
  real function UserUPlus(grid, y)
!----------------------------------------------------------------------!
! Prescribes initial velocity profile for channel flow at ReTau=180    !
! This function gives centerline velocity: Uc = 18.0412006             !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  real :: y
!======================================================================!
!   Wall is at y=0, and centerline at y=1                              !
!----------------------------------------------------------------------!

  UserUPlus =                                                       &
   -4.555442e+003  * y**9  +  1.902458e+004  * y**8                 &
   -3.108515e+004  * y**7  +  2.331028e+004  * y**6                 &
   -4.723039e+003  * y**5  -  4.867699e+003  * y**4                 &
   +3.999634e+003  * y**3  -  1.297551e+003  * y**2                 &
   +2.125647e+002  * y**1  -  1.353124e-001

  end function UserUPlus
