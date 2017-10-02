!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
MODULE par_mod

  use allp_mod

  implicit none

  integer :: Nsub, subNC(MAXPRO), NBBs(0:MAXPRO), NBBe(0:MAXPRO),   &
             this, NPro
  integer,allocatable :: proces(:), BuSeIn(:), BuReIn(:), &
                         BufInd(:), BufPos(:) 

end MODULE
