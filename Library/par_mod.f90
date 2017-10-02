!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
MODULE par_mod

  USE allp_mod

  IMPLICIT NONE

  INTEGER :: Nsub, subNC(MAXPRO), NBBs(0:MAXPRO), NBBe(0:MAXPRO),   &
             this, NPro
  INTEGER,ALLOCATABLE :: proces(:), BuSeIn(:), BuReIn(:), &
                         BufInd(:), BufPos(:) 

END MODULE
