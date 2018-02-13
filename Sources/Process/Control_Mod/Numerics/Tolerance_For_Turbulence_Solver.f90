!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Turbulence_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_FOR_TURBULENCE_SOLVER', 1.0e-3,  &
                                   val, verbose)

  end subroutine
