!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE div_mod

  implicit none

  integer,allocatable :: ix(:), iy(:), iz(:), iin(:)
  real,allocatable    :: criter(:) 

!------------!
! Parameters !
!------------!
  integer :: ALGOR, COORDINATE, INERTIAL

end MODULE
