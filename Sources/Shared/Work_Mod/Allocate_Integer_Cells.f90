!==============================================================================!
  subroutine Work_Mod_Allocate_Integer_Cells(grid)
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc, nb
!==============================================================================!

  ! Get number of cells and boundary cells
  nc = grid % n_cells
  nb = grid % n_bnd_cells

  end subroutine
