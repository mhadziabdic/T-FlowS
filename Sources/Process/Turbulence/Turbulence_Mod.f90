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

  end module 
