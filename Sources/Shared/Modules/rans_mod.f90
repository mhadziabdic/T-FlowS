!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   RANS models  variable         !   Delft University of Technology   !
!   definitions for the processor !   Section Heat Transfer             !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

module rans_mod

  use allt_mod, only: Unknown

  implicit none

!----- Turbulence models variables
  type(Unknown) :: KIN
  type(Unknown) :: EPS
  type(Unknown) :: V_2
  type(Unknown) :: F22
  type(Unknown) :: VIS

  type(Unknown) :: uu
  type(Unknown) :: vv
  type(Unknown) :: ww
  type(Unknown) :: uv
  type(Unknown) :: uw
  type(Unknown) :: vw
 
!----- Constants for the k-eps model:
  real :: Ce1, Ce2, Ce3, Cmu, Cmu25, Cmu75, kappa, Elog, Zo
 
!----- Constants for the k-eps-v2f model:
  real :: CmuD, Cl, Ct, alpha, Cni, cf1, cf2, cf3, Cf_1, Cf_2
  real :: Lim
  real :: g1, g1_star, g2, g3, g3_star, g4, g5 

!----- Constants for the Spalart-Allmaras model:
  real :: Cb1, Cb2, SIGMAv, Cw1, Cw2, Cw3, Cvis1

!----- Vorticity
  real,allocatable :: Vort(:), VortMean(:)

!----- Turbulent viscosity
  real,allocatable :: VISt(:), CmuS(:)
 
!----- Turbulent conductivity
  real,allocatable :: CONt(:)
 
!----- Lenght and Time Scales
  real,allocatable :: Lsc(:)
  real,allocatable :: Tsc(:)   

!----- Production of turbulent kinetic energy
  real,allocatable :: Pk(:)
!----- Buoyancy production
  real,allocatable :: Gbuoy(:)
  real,allocatable :: buoyBeta(:)
  real,allocatable :: Pbuoy(:)
 
!----- Non-dimensional distance
  real,allocatable :: Ynd(:)
 
!----- Friction velocity
  real,allocatable :: Uf(:)
  real,allocatable :: Ufmean(:)

!----- Gravity
  real :: grav_x, grav_y, grav_z

!----- Wall viscosity (wall function approuch)
  real,allocatable :: VISwall(:)
  real,allocatable :: CONwall(:)

!  real,allocatable :: AA(:)
!  real,allocatable :: EE(:)
  real,allocatable :: Fs(:)
!  real,allocatable :: Feps(:)

  real,allocatable :: nn1(:)
  real,allocatable :: nn2(:)
  real,allocatable :: nn3(:)

  real,allocatable :: Bud1(:)
  real,allocatable :: Bud2(:)
  real,allocatable :: Bud3(:)
  real,allocatable :: Bud4(:)
  real,allocatable :: Bud5(:)
  real,allocatable :: Bud6(:)
  real,allocatable :: Bud7(:)
  real,allocatable :: Bud8(:)
  real,allocatable :: Bud9(:)
 
  real,allocatable :: uu_star(:)
  real,allocatable :: vv_star(:)
  real,allocatable :: ww_star(:)
  real,allocatable :: uv_star(:)
  real,allocatable :: uw_star(:)
  real,allocatable :: vw_star(:)

end module 
