!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!   Domain as the one used in "Generator"                                      !
!------------------------------------------------------------------------------!
  use Line_Mod
  use Block_Mod
  use Range_Mod                ! to store boundary conditions or materials
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Domain type   !
  !-----------------!
  type Domain_Type

    type(Block_Type),              allocatable :: blocks(:)
    type(Line_Type),               allocatable :: lines(:)
    type(Range_Type),              allocatable :: ranges(:)

  end type Domain_Type

  ! If defined like this, one can easily think of multiple domains
  type(Domain_Type) :: dom

  end module Domain_Mod
