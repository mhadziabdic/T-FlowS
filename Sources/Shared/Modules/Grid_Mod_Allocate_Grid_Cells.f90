!==============================================================================!
  subroutine Allocate_Grid_Cells(grid, nb, ni)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nb        ! number of cells on the bounday
  integer         :: ni        ! number of cells inside
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  allocate(grid % cells(-nb:ni))

  ! Initialize nodes and cells
  do c = -nb, ni
    grid % cells(c) % n_nodes = 0 
    grid % cells(c) % n       = 0 
    grid % cells(c) % c       = 0 
  end do

  ! Initialize coordinates
  grid % cells(-nb:ni) % x = 0 
  grid % cells(-nb:ni) % y = 0 
  grid % cells(-nb:ni) % z = 0 

  end subroutine
