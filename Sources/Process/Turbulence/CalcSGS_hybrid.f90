!======================================================================!
  subroutine CalcSGS_hybrid()
!----------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for LES.           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!----------------------------------------------------------------------!
! near(c) is the number of corresponding cell on the nearest wall.
! In case that, in parallel executions, the subdomain does not have 
! any nearwall cells, the near(c) is zero.
! near(c) is calculated in NearWallCells.f90, only ones in the beginig
! of a simulation.
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, s, c1, c2 
  real    :: Nx, Ny, Nz
  real    :: Cs, R
  real    :: Stot, lf, UtauL, Uff 
  real    :: Utot, Unor, Utan, Apow, Bpow, nu, dely, yPlus 
  real    :: fun
!======================================================================!

  
!===================!
!                   !
!     SGS terms     !
!                   !
!===================!

  do c=1,NC
    lf =  (volume(c)**ONE_THIRD)    
    VISt_sgs(c) = DENc(material(c))    &
                  * (lf*lf)              &          ! delta^2 
                  * Cdyn(c)              &          ! Cdynamic   
                  * Shear(c)      
  end do

  call Exchng(VISt_sgs)

  end subroutine
