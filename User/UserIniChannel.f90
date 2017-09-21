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
!--------------------------------[CVS]---------------------------------!
!  $Id: UserIniChannel.f90,v 1.1 2017/08/31 22:42:35 mhadziabdic Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/User/UserIniChannel.f90,v $  
!======================================================================!

    do c=1,NC
      yPlus = WallDs(c)*590
      U% n(c) = (yPlus + 1.0/0.41 * (LOG())
    end do
 


  END SUBROUTINE UserIniChannel
