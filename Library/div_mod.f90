!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE div_mod

  IMPLICIT NONE

  INTEGER,ALLOCATABLE :: ix(:), iy(:), iz(:), iin(:)
  REAL,ALLOCATABLE    :: criter(:) 

!------------!
! Parameters !
!------------!
  INTEGER :: ALGOR, COORDINATE, INERTIAL

END MODULE
