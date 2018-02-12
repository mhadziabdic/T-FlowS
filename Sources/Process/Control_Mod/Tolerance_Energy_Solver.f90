!==============================================================================!
  subroutine Control_Mod_Tolerance_Energy_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_ENERGY_SOLVER', 1.0e-3,  &
                                   val, verbose)

  end subroutine
