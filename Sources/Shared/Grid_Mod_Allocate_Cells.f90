!==============================================================================!
  subroutine Grid_Mod_Allocate_Cells(grid, nb, ni)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nb        ! number of cells on the bounday
  integer         :: ni        ! number of cells inside
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Allocate cell center coordinates and initialize to zero
  allocate(grid % xc(-nb:ni));  grid % xc = 0.0
  allocate(grid % yc(-nb:ni));  grid % yc = 0.0
  allocate(grid % zc(-nb:ni));  grid % zc = 0.0

  ! Cells' nodes and neigboring cells
  allocate(grid % cells_n( 8, -nb:ni));    grid % cells_n       = 0
  allocate(grid % cells_c(24, -nb:ni));    grid % cells_c       = 0

  ! Number of nodes at each cell (determines cell's shape really)
  allocate(grid % cells_n_nodes(-nb:ni));  grid % cells_n_nodes = 0

  end subroutine
