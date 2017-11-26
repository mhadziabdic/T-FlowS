!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                                  !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
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
  allocate (delta (-grid % n_bnd_cells:grid % n_cells));  delta  = 0.0  
  allocate (WallDs(-grid % n_bnd_cells:grid % n_cells));  WallDs = 0.0       
  allocate (fF(grid % n_faces));  fF = 0.0  

  ! Variables defined in Solvers_Mod
  call Matrix_Mod_Allocate(grid, D)

  ! Variables defined in pro_mod.h90:
  call Matrix_Mod_Allocate(grid, A)
  allocate (b(grid % n_cells));  b=0

  ! Working arrays
  call Work_Mod_Allocate_Real_Cells(grid, 15)
  call Work_Mod_Allocate_Real_Faces(grid,  1)

  allocate (Scoef(grid % n_faces)); Scoef=0.

  allocate (xp(grid % n_materials));    xp   =0.0
  allocate (yp(grid % n_materials));    yp   =0.0
  allocate (zp(grid % n_materials));    zp   =0.0

  ! Variables defined in par_mod.h90:
  allocate (BufInd(-grid % n_bnd_cells:-1)); BufInd=0

  end subroutine
