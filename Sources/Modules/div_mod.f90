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

  integer, parameter :: COORDINATE = 1
  integer, parameter :: INERTIAL   = 2

  integer              :: division_algorithm
  integer, allocatable :: ix(:), iy(:), iz(:), iin(:)
  real, allocatable    :: criter(:) 

end MODULE
