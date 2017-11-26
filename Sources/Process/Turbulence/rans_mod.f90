!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   RANS models  variable         !   Delft University of Technology   !
!   definitions for the processor !   Section Heat Transfer             !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

  module rans_mod

  use Var_Mod

  implicit none

  ! Turbulence models variables
  type(Var_Type) :: kin
  type(Var_Type) :: eps
  type(Var_Type) :: v_2
  type(Var_Type) :: f22
  type(Var_Type) :: vis

  ! Reynolds stresses
  type(Var_Type) :: uu
  type(Var_Type) :: vv
  type(Var_Type) :: ww
  type(Var_Type) :: uv
  type(Var_Type) :: uw
  type(Var_Type) :: vw
 
  ! Temperature fluctuations
  type(Var_Type) :: tt
  type(Var_Type) :: ut
  type(Var_Type) :: vt
  type(Var_Type) :: wt
 
  ! Constants for the k-eps model:
  real :: Ce1, Ce2, Ce3, Cmu, Cmu25, Cmu75, kappa, Elog, Zo
 
  ! Constants for the k-eps-v2f model:
  real :: CmuD, Cl, Ct, alpha, Cni, cf1, cf2, cf3, Cf_1, Cf_2
  real :: Lim
  real :: g1, g1_star, g2, g3, g3_star, g4, g5 

  ! Constants for the Spalart-Allmaras model:
  real :: Cb1, Cb2, SIGMAv, Cw1, Cw2, Cw3, Cvis1

  ! Total dissipation in HJ model
  real,allocatable :: eps_tot(:)

  ! Vorticity
  real,allocatable :: Vort(:), VortMean(:)

  ! Turbulent viscosity
  real,allocatable :: VISt(:), CmuS(:)
 
  ! Turbulent conductivity
  real,allocatable :: CONt(:)
 
  ! Lenght and Time Scales
  real,allocatable :: Lsc(:)
  real,allocatable :: Tsc(:)   

  ! Production of turbulent kinetic energy
  real,allocatable :: Pk(:)
  ! Buoyancy production
  real,allocatable :: Gbuoy(:)
  real,allocatable :: buoyBeta(:)
  real,allocatable :: Pbuoy(:)
 
  ! Non-dimensional distance
  real,allocatable :: Ynd(:)
 
  ! Friction velocity
  real,allocatable :: Uf(:)
  real,allocatable :: Ufmean(:)

  ! Gravity
  real :: grav_x, grav_y, grav_z

  ! Wall viscosity (wall function approuch)
  real,allocatable :: VISwall(:)
  real,allocatable :: CONwall(:)

  real,allocatable :: Fs(:)

  ! These are some working variables for RSM models
  real,allocatable :: VAR1x(:),   VAR1y(:),   VAR1z(:)
  real,allocatable :: VAR2x(:),   VAR2y(:),   VAR2z(:)
  real,allocatable :: VAR3x(:),   VAR3y(:),   VAR3z(:)
  real,allocatable :: VAR4x(:),   VAR4y(:),   VAR4z(:)
  real,allocatable :: VAR5x(:),   VAR5y(:),   VAR5z(:)
  real,allocatable :: VAR6x(:),   VAR6y(:),   VAR6z(:)
  real,allocatable :: VAR7x(:),   VAR7y(:),   VAR7z(:)
  real,allocatable :: VAR8x(:),   VAR8y(:),   VAR8z(:)
  real,allocatable :: VAR9x(:),   VAR9y(:),   VAR9z(:)
  real,allocatable :: VAR10x(:),  VAR10y(:),  VAR10z(:)
  real,allocatable :: VAR11x(:),  VAR11y(:),  VAR11z(:)
  real,allocatable :: VAR12x(:),  VAR12y(:),  VAR12z(:)
  real,allocatable :: PHI1x(:),   PHI1y(:),   PHI1z(:)
  real,allocatable :: PHI2x(:),   PHI2y(:),   PHI2z(:)
  real,allocatable :: PHI3x(:),   PHI3y(:),   PHI3z(:)
  real,allocatable :: PHI4x(:),   PHI4y(:),   PHI4z(:)
  real,allocatable :: PHI5x(:),   PHI5y(:),   PHI5z(:)
  real,allocatable :: PHI6x(:),   PHI6y(:),   PHI6z(:)
  real,allocatable :: PHI7x(:),   PHI7y(:),   PHI7z(:)
  real,allocatable :: PHI8x(:),   PHI8y(:),   PHI8z(:)
  real,allocatable :: PHI9x(:),   PHI9y(:),   PHI9z(:)
  real,allocatable :: PHI10x(:),  PHI10y(:),  PHI10z(:)
  real,allocatable :: PHI11x(:),  PHI11y(:),  PHI11z(:)

  end module 
