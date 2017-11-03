!==============================================================================!
  subroutine Var_Mod_Allocate_New(phi, n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: n_bnd_cells
  integer        :: n_cells
!------------------------------------------------------------------------------!

  ! Values in the new (n) time step
  allocate (phi % n (-n_bnd_cells: n_cells));   phi % n  = 0.

  end subroutine
