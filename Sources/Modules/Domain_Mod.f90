!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!   Domain as the one used in "Generator"                                      !
!------------------------------------------------------------------------------!
  use Line_Mod
  use Block_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Domain type   !
  !-----------------!
  type Domain_Type

    type(Block_Type), allocatable :: blocks(:)
    type(Line_Type),  allocatable :: lines(:)

  end type Domain_Type

  ! If defined like this, one can easily think of multiple domains
  type(Domain_Type) :: dom

  end module Domain_Mod
