!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
module par_mod

  use allp_mod

  implicit none

  integer :: n_sub      ! Number of subdivisions.  Used in Divide
  integer :: this_proc  ! Processor i.d.
  integer :: n_proc     ! Number of processors.  Mostly in Process

  integer, allocatable :: subNC(:), NBBs(:), NBBe(:)

  integer, allocatable :: proces(:), BuSeIn(:), BuReIn(:), &
                          BufInd(:), BufPos(:) 

end module
