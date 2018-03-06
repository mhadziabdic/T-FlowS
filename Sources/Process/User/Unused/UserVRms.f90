!======================================================================!
  real function UserVRms(grid, y)
!----------------------------------------------------------------------!
!   Prescribes mean spanwise velocity fluctuations for the channel.    !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  real :: y
!======================================================================!
!   Wall is at y=0, and centerline at y=1                              !
!----------------------------------------------------------------------!

  UserVRms =                                                        &
   +2.432998e+003  * y**9  -  1.166312e+004  * y**8                 &
   +2.368065e+004  * y**7  -  2.653575e+004  * y**6                 &
   +1.793421e+004  * y**5  -  7.514617e+003  * y**4                 &
   +1.941100e+003  * y**3  -  3.012721e+002  * y**2                 &
   +2.639077e+001  * y**1  +  2.737377e-002

  end function UserVRms