!==============================================================================!
  module Material_Mod
!------------------------------------------------------------------------------!
!   This is used to store materials in "Generator"                             !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Material type   !
  !-------------------!
  type Material_Type

    ! Stores the name of the material
    character(len=80) :: name

  end type Material_Type

  end module Material_Mod
