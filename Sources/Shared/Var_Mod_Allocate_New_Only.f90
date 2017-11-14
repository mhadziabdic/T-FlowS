!==============================================================================!
  subroutine Var_Mod_Allocate_New_Only(name_phi, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: name_phi
  type(Var_Type)   :: phi
  type(Grid_Type)  :: grid
!==============================================================================!

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n (-grid % n_bnd_cells : grid % n_cells));   phi % n = 0.

  end subroutine
