!======================================================================!
  subroutine UserIniChannel
!----------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use par_mod
  use les_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c1, c2, s, m, step
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
  real    :: Stot
  real    :: AreaCross, AreaLowWall, AreaHighWall
  real    :: Tbulk, TlowWall, ThighWall
!======================================================================!

    do c=1,NC
      yPlus = WallDs(c)*590
      U% n(c) = (yPlus + 1.0/0.41 * (LOG())
    end do
 


  end subroutine UserIniChannel
