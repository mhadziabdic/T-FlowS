!==============================================================================!
  subroutine Grid_Mod_Allocate_Cells(grid, nb, nc)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nb        ! number of cells on the bounday
  integer         :: nc        ! number of cells inside
!==============================================================================!

  ! Store number of cells and boundary cells
  grid % n_cells          = nc
  grid % n_bnd_cells = nb

  ! Allocate cell center coordinates and initialize to zero
  allocate(grid % xc(-nb:nc));  grid % xc = 0.0
  allocate(grid % yc(-nb:nc));  grid % yc = 0.0
  allocate(grid % zc(-nb:nc));  grid % zc = 0.0

  ! Memory for cells' volumes
  allocate(grid % vol(-nb:nc));  grid % vol = 0.0

  ! Cells' nodes and neigboring cells
  allocate(grid % cells_n( 8, -nb:nc));    grid % cells_n       = 0
  allocate(grid % cells_c(24, -nb:nc));    grid % cells_c       = 0

  ! Number of nodes at each cell (determines cell's shape really)
  allocate(grid % cells_n_nodes(-nb:nc));  grid % cells_n_nodes = 0

  end subroutine
