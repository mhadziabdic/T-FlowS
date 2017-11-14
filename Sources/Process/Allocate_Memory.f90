!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                                  !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use Grid_Mod
  use Matrix_Mod
  use Solvers_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Variables defined in all.h90:
  allocate (delta (-grid % n_boundary_cells:grid % n_cells));  delta  = 0.0  
  allocate (WallDs(-grid % n_boundary_cells:grid % n_cells));  WallDs = 0.0       
  allocate (f (grid % n_faces));  f  = 0.0  
  allocate (fF(grid % n_faces));  fF = 0.0  

  ! Variables defined in sol.h90:
  call Solvers_Mod_Allocate_Vectors(grid % n_boundary_cells, grid % n_cells)
  call Matrix_Mod_Allocate(D, grid % n_boundary_cells,  &
                              grid % n_cells,           &
                              grid % n_faces)

  ! Variables defined in pro_mod.h90:
  call Matrix_Mod_Allocate(A, grid % n_boundary_cells,  &
                              grid % n_cells,           &
                              grid % n_faces)
  allocate (b(grid % n_cells));  b=0

  allocate (Scoef(grid % n_faces)); Scoef=0.

  allocate (xp(grid % n_materials));    xp   =0.0
  allocate (yp(grid % n_materials));    yp   =0.0
  allocate (zp(grid % n_materials));    zp   =0.0
  allocate (AreaX(grid % n_materials)); AreaX=0.0
  allocate (AreaY(grid % n_materials)); AreaY=0.0
  allocate (AreaZ(grid % n_materials)); AreaZ=0.0

  ! Variables defined in par_mod.h90:
  allocate (BufInd(-grid % n_boundary_cells:-1)); BufInd=0

  end subroutine
