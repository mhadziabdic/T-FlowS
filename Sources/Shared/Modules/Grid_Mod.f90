!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids module is used throughout all programs                               !
!   (that means in "Generate", "Divide", "Neu2TflowS", "Process".              !
!------------------------------------------------------------------------------!
  use Material_Mod
  use Boundary_Condition_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Grid type   !
  !               !
  !---------------!
  type Grid_Type

    integer :: n_materials            ! number of materials
    integer :: n_boundary_conditions  ! number of boundary conditions

    !-------------------------!
    !  Cell-based variables   !
    !-------------------------!

    ! Cell center coordinates
    real, allocatable :: xc(:), yc(:), zc(:)  
    
    ! Cells' nodes and neigboring cells
    integer, allocatable :: cells_n(:,:)      
    integer, allocatable :: cells_c(:,:)

    ! Number of nodes at each cell (determines cell's shape really)
    integer, allocatable :: cells_n_nodes(:)

    !-------------------------!
    !  Face-based variables   !
    !-------------------------!

    ! Number of nodes at each face (determines face's shape really)
    integer, allocatable :: faces_n_nodes(:)

    ! Faces' nodes and neigboring cells
    integer, allocatable :: faces_n(:,:)
    integer, allocatable :: faces_c(:,:)

    !-------------------------!
    !  Node-based variables   !
    !-------------------------!

    ! Node coordinates
    real, allocatable :: xn(:), yn(:), zn(:)
    
    type(Material_Type),           allocatable :: materials(:)
    type(Boundary_Condition_Type), allocatable :: boundary_conditions(:)

    !  Maximum number of cells, boundary cells and faces
    ! (Used for tentative memory allocation in Generator)
    integer :: max_n_nodes
    integer :: max_n_boundary_cells
    integer :: max_n_faces

  end type Grid_Type

  ! If defined like this, one can easily think of multiple grids
  ! type(Grid_Type) :: grid

  contains
 
  include 'Grid_Mod_Sort_Cells_By_Index.f90'
  include 'Grid_Mod_Sort_Faces_By_Index.f90'
  include 'Grid_Mod_Allocate_Cells.f90'
  include 'Grid_Mod_Allocate_Faces.f90'
  include 'Grid_Mod_Allocate_Nodes.f90'

  end module Grid_Mod
