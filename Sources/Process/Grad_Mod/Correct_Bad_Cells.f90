!==============================================================================!
  subroutine Grad_Mod_Correct_Bad_Cells(grid, phii)
!------------------------------------------------------------------------------!
!   Corrects the pressure gradients in the cells where they cannot             !
!   be computed, the so called "bad" cells.                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phii(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
!==============================================================================!

  do c = 1, grid % n_cells
    if(bad_cells(c)) then
      phii(c) = 0.0
    end if
  end do 

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
     
    if(c2 > 0 .or.  &
       c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then
      if(bad_cells(c1)) phii(c1) = phii(c1) + 0.5*phii(c2) 
      if(bad_cells(c2)) phii(c2) = phii(c2) + 0.5*phii(c1) 
    end if
  end do

  call Comm_Mod_Exchange(grid, phii)

  end subroutine