!==============================================================================!
  subroutine Grid_Mod_Allocate_Faces(grid, nf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nf        ! number of faces in the grid   
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Number of nodes at each face (determines face's shape really)
  allocate(grid % faces_n_nodes(nf));  grid % faces_n_nodes = 0

  ! Faces' nodes and neigboring cells
  allocate(grid % faces_n(4, nf));  grid % faces_n = 0
  allocate(grid % faces_c(2, nf));  grid % faces_c = 0

  end subroutine
