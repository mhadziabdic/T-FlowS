!======================================================================!
  SUBROUTINE UserIniChannel
!----------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c1, c2, s, m, step
  REAL    :: xc1, yc1, zc1, xc2, yc2, zc2
  REAL    :: Stot
  REAL    :: AreaCross, AreaLowWall, AreaHighWall
  REAL    :: Tbulk, TlowWall, ThighWall
!======================================================================!

    do c=1,NC
      yPlus = WallDs(c)*590
      U% n(c) = (yPlus + 1.0/0.41 * (LOG())
    end do
 


  END SUBROUTINE UserIniChannel
