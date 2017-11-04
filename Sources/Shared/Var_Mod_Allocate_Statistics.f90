!==============================================================================!
  subroutine Var_Mod_Allocate_Statistics(phi, n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
!   This is to allocate additional values for statistics.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: n_bnd_cells
  integer        :: n_cells
!------------------------------------------------------------------------------!

  ! Terms for statistics
  allocate (phi % mean(-n_bnd_cells: n_cells));  phi % mean = 0.
  allocate (phi % fluc(-n_bnd_cells: n_cells));  phi % fluc = 0.
  allocate (phi % filt(-n_bnd_cells: n_cells));  phi % filt = 0.

  end subroutine
