!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                                  !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use par_mod
  use Grid_Mod
  use Work_Mod
  use Matrix_Mod
  use Solvers_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Variables defined in all.h90:
  allocate (fw(grid % n_faces));  fw = 0.0  

  ! Variables defined in Solvers_Mod
  call Matrix_Mod_Allocate(grid, D)

  ! Variables defined in Flow_Mod.h90:
  call Matrix_Mod_Allocate(grid, A)
  allocate (b(grid % n_cells));  b=0

  ! Working arrays
  call Work_Mod_Allocate_Real_Cells(grid, 20)
  call Work_Mod_Allocate_Real_Faces(grid,  1)
  call Work_Mod_Allocate_Real_Nodes(grid,  1)

  ! This array should be handled in a more elegant way
  allocate (f_coef(grid % n_faces)); f_coef=0.

  ! Variables defined in par_mod.h90:
  allocate (BufInd(-grid % n_bnd_cells:-1)); BufInd=0

  ! Variables actively used in Cgns_Mod
 
  end subroutine
