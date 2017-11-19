!==============================================================================!
  module Boundary_Condition_Mod
!------------------------------------------------------------------------------!
!   This is used to store boundary conditions or materials in "Generator"      !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------------------!
  !   Boundary_Conditions type   !
  !------------------------------!
  type Boundary_Condition_Type

    character(len=80) :: name

  end type

  end module
