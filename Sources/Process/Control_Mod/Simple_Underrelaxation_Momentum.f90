!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_Momentum(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SIMPLE_UNDERRELAXATION_MOMENTUM', 0.7,  &
                                   val, verbose)

  end subroutine
