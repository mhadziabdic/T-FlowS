!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module div_mod

  implicit none

  ! Division algorithm
  integer            :: division_algorithm
  integer, parameter :: COORDINATE = 1
  integer, parameter :: INERTIAL   = 2

  ! Number of sub-divisions
  integer              :: n_sub

  ! Buffer send index and buffer receive index.  
  ! Used for plotting dcomposed grids with links.
  integer, allocatable :: SubNC(:), BuSeIn(:), BuReIn(:), BufPos(:)

  ! Processor i.d.
  integer, allocatable :: proces(:)  

  ! Axes of innertial
  integer, allocatable :: ix(:), iy(:), iz(:), iin(:)

  ! Division criterion (I believe)
  real, allocatable    :: criter(:) 

end module
