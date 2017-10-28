!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids defining blocks used in "Generator"                                  !
!------------------------------------------------------------------------------!
  use Point_Mod
  use Cell_Mod
  use Material_Mod
  use Boundary_Condition_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Grid type   !
  !---------------!
  type Grid_Type

    ! Maximum number of cells, boundary cells and faces
    !  (They are used for tentative memory allocation)
    integer :: max_n_nodes
    integer :: max_n_boundary_cells
    integer :: max_n_faces

    type(Point_Type),              allocatable :: nodes(:)
    type(Cell_Type),               allocatable :: cells(:)  

    type(Material_Type),           allocatable :: materials(:)
    type(Boundary_Condition_Type), allocatable :: boundary_conditions(:)

  end type Grid_Type

  ! If defined like this, one can easily think of multiple grids
  type(Grid_Type) :: grid

  contains
 
  include 'Grid_Mod_Sort_Cells_By_Index.f90'
  include 'Grid_Mod_Allocate_Grid_Cells.f90'
  include 'Grid_Mod_Allocate_Grid_Nodes.f90'

  end module Grid_Mod
