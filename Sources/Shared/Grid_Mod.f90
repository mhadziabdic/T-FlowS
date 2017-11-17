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

    ! Number of ...
    integer :: n_nodes                ! ... nodes
    integer :: n_cells                ! ... cells
    integer :: n_faces                ! ... faces
    integer :: n_bnd_cells            ! ... boundary cells
    integer :: n_per_faces            ! ... periodic faces (shadows)
    integer :: n_materials            ! ... materials
    integer :: n_boundary_conditions  ! ... boundary conditions
    integer :: n_copy                 ! ... copy cells and faces
    integer :: n_sh                   ! ... shadow faces           

    !-------------------------!
    !  Cell-based variables   !
    !-------------------------!

    ! Cell center coordinates
    real, allocatable :: xc(:), yc(:), zc(:)  

    ! Cell volumes
    real, allocatable :: vol(:)
    
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

    ! Face surface areas (si) and distances between cells (di)
    real, allocatable :: sx(:), sy(:), sz(:)
    real, allocatable :: dx(:), dy(:), dz(:)

    ! Face coordinates 
    real, allocatable :: xf(:), yf(:), zf(:)
    
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
    integer :: max_n_bnd_cells
    integer :: max_n_faces

  end type

  ! If defined like this, one can easily think of multiple grids
  ! type(Grid_Type) :: grid

  contains
 
  include 'Grid_Mod/Sort_Cells_By_Index.f90'
  include 'Grid_Mod/Sort_Faces_By_Index.f90'
  include 'Grid_Mod/Allocate_Cells.f90'
  include 'Grid_Mod/Allocate_Faces.f90'
  include 'Grid_Mod/Allocate_Nodes.f90'

  end module
