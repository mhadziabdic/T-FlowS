!==============================================================================!
  subroutine Work_Mod_Allocate_Integer_Cells(grid, n)
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc, nb
!==============================================================================!

  ! Get number of cells and boundary cells
  nc = grid % n_cells
  nb = grid % n_bnd_cells

  ! Store the pointer to the grid
  pnt_grid => grid

  ! Allocate requested memory
  allocate(i_cell_01(-nb:nc));   i_cell_01 = 0;   if(n ==  1) return
  allocate(i_cell_02(-nb:nc));   i_cell_02 = 0;   if(n ==  2) return
  allocate(i_cell_03(-nb:nc));   i_cell_03 = 0;   if(n ==  3) return
  allocate(i_cell_04(-nb:nc));   i_cell_04 = 0;   if(n ==  4) return
  allocate(i_cell_05(-nb:nc));   i_cell_05 = 0;   if(n ==  5) return
  allocate(i_cell_06(-nb:nc));   i_cell_06 = 0;   if(n ==  6) return
  allocate(i_cell_07(-nb:nc));   i_cell_07 = 0;   if(n ==  7) return
  allocate(i_cell_08(-nb:nc));   i_cell_08 = 0;   if(n ==  8) return
  allocate(i_cell_09(-nb:nc));   i_cell_09 = 0;   if(n ==  9) return
  allocate(i_cell_10(-nb:nc));   i_cell_10 = 0;   if(n == 10) return
  allocate(i_cell_11(-nb:nc));   i_cell_11 = 0;   if(n == 11) return
  allocate(i_cell_12(-nb:nc));   i_cell_12 = 0;   if(n == 12) return
  allocate(i_cell_13(-nb:nc));   i_cell_13 = 0;   if(n == 13) return
  allocate(i_cell_14(-nb:nc));   i_cell_14 = 0;   if(n == 14) return
  allocate(i_cell_15(-nb:nc));   i_cell_15 = 0;   if(n == 15) return
  allocate(i_cell_16(-nb:nc));   i_cell_16 = 0;   if(n == 16) return
  allocate(i_cell_17(-nb:nc));   i_cell_17 = 0;   if(n == 17) return
  allocate(i_cell_18(-nb:nc));   i_cell_18 = 0;   if(n == 18) return
  allocate(i_cell_19(-nb:nc));   i_cell_19 = 0;   if(n == 19) return
  allocate(i_cell_20(-nb:nc));   i_cell_20 = 0;   if(n == 20) return
  allocate(i_cell_21(-nb:nc));   i_cell_21 = 0;   if(n == 21) return
  allocate(i_cell_22(-nb:nc));   i_cell_22 = 0;   if(n == 22) return
  allocate(i_cell_23(-nb:nc));   i_cell_23 = 0;   if(n == 23) return
  allocate(i_cell_24(-nb:nc));   i_cell_24 = 0;   if(n == 24) return
  allocate(i_cell_25(-nb:nc));   i_cell_25 = 0;   if(n == 25) return
  allocate(i_cell_26(-nb:nc));   i_cell_26 = 0;   if(n == 26) return
  allocate(i_cell_27(-nb:nc));   i_cell_27 = 0;   if(n == 27) return
  allocate(i_cell_28(-nb:nc));   i_cell_28 = 0;   if(n == 28) return
  allocate(i_cell_29(-nb:nc));   i_cell_29 = 0;   if(n == 29) return
  allocate(i_cell_30(-nb:nc));   i_cell_30 = 0;   if(n == 30) return

  end subroutine
