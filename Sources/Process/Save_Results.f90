!==============================================================================!
  subroutine Save_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   If T-Flows was compiled with CGNS_ADF5=yes or CGNS_HDF5=yes                !
!   Then use Save_Cgns_Results to save fields                                  !
!   Otherwise Save_Vtu_Results is used                                         !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!==============================================================================!

  if (COMPILED_WITH_CGNS) then
    call Save_Cgns_Results(grid, name_save)
  else
    call Save_Vtu_Results (grid, name_save)
  end if

  end subroutine