!==============================================================================!
  module Turbulence_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all turbulence modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Varibale holing the turbulence model
  integer :: turbulence_model

  ! Parameters describing turbulence model choice
  integer, parameter :: K_EPS                 = 30011
  integer, parameter :: K_EPS_V2              = 30013
  integer, parameter :: K_EPS_ZETA_F          = 30029
  integer, parameter :: HYBRID_K_EPS_ZETA_F   = 30047
  integer, parameter :: HYBRID_PITM           = 30059
  integer, parameter :: LES                   = 30071
  integer, parameter :: DNS                   = 30089
  integer, parameter :: DES_SPALART           = 30091
  integer, parameter :: SPALART_ALLMARAS      = 30097  
  integer, parameter :: HANJALIC_JAKIRLIC     = 30103
  integer, parameter :: REYNOLDS_STRESS_MODEL = 30109

  ! Varibale holing turbulence model variant
  integer :: turbulence_model_variant

  ! Turbulence model variants
  integer, parameter :: NONE        = 30113
  integer, parameter :: HYBRID      = 30119
  integer, parameter :: PURE        = 30133
  integer, parameter :: URANS       = 30137
  integer, parameter :: LOW_RE      = 30139
  integer, parameter :: HIGH_RE     = 30161
  integer, parameter :: WALE        = 30169
  integer, parameter :: DYNAMIC     = 30181
  integer, parameter :: SMAGORINSKY = 30187

  ! Reynolds stresses
  type(Var_Type) :: uu
  type(Var_Type) :: vv
  type(Var_Type) :: ww
  type(Var_Type) :: uv
  type(Var_Type) :: uw
  type(Var_Type) :: vw
 
  ! Turbulent heat fluxes
  type(Var_Type) :: tt
  type(Var_Type) :: ut
  type(Var_Type) :: vt
  type(Var_Type) :: wt
 
  ! Tripple moments
  type(Var_Type) :: uuu, uuv, uuw
  type(Var_Type) :: vvu, vvv, vvw
  type(Var_Type) :: wwu, wwv, www
  type(Var_Type) :: uwu, uwv, uww

  ! Shear and wall stress are used in a number of turbulence models
  real, allocatable :: shear(:)
  real, allocatable :: tau_wall(:)

  ! Wall viscosity (wall function approuch)
  real, allocatable :: vis_wall(:)
  real, allocatable :: con_wall(:)

  ! Non-dimensional distance
  real, allocatable :: y_plus(:)
 
  ! Friction velocity and its time-average
  real,allocatable :: u_tau(:)
  real,allocatable :: u_tau_mean(:)

  end module 
