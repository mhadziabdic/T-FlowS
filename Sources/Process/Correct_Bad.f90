!==============================================================================!
  subroutine Correct_Bad(grid, phii)
!------------------------------------------------------------------------------!
!   Corrects the pressure gradients in the cells where they cannot             !
!   be computed, the so called "bad" cells.                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
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
    if(BadForG(c)) then
      phii(c) = 0.0
    end if
  end do 

  do s = 1, grid % n_faces
    c1 = SideC(1,s)
    c2 = SideC(2,s)
     
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then
      if(BadForG(c1)) phii(c1) = phii(c1) + 0.5*phii(c2) 
      if(BadForG(c2)) phii(c2) = phii(c2) + 0.5*phii(c1) 
    end if
  end do

  call Exchng(phii)

  end subroutine
