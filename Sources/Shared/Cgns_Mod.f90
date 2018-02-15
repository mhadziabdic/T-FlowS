!==============================================================================!
  module Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!==============================================================================!

  !----------!
  !   File   ! -> contains base
  !----------!
    !
    !----------!
    !   Base   ! -> contains blocks
    !----------!
      !
      !------------!
      !   Blocks   ! -> contains coordinates, sections, bnd_conds, solution
      !------------!
        !
        !-----------------!
        !   Coordinates   !
        !-----------------!
        !
        !----------------------!
        !   Element sections   !
        !----------------------!
        !
        !-------------------------!
        !   Boundary conditions   !
        !-------------------------!
        !
        !----------------------!
        !   Solution section   !
        !----------------------!
          !
          !------------------!
          !   Field section  !
          !------------------!

! TO DO : Pack monitoring points in ConvergenceHistory_t

  ! File
  integer           :: file_id
  character(len=80) :: file_name
  integer           :: file_mode
  logical           :: verbose = .true.

  ! Field section
  type Cgns_Field_Type
    character(len=80) :: name
    integer           :: field_type ! RealDouble & LongInt?
  end type

  ! Solution section
  type Cgns_Solution_Type
    character(len=80)                  :: name
    integer                            :: sol_type
    integer                            :: n_fields
    type(Cgns_Field_Type), allocatable :: field(:)
  end type

  ! Element section
  type Cgns_Section_Type
    character(len=80)    :: name
    integer, allocatable :: cell_type
    integer              :: first_cell
    integer              :: last_cell
    integer              :: parent_flag
  end type

  ! Boundary conditions   ! -> it is similar to Bnd_Cond in ../Share :-(
  type Cgns_Bnd_Cond_Type
    character(len=80)    :: name
    integer, allocatable :: color
    integer, allocatable :: n_nodes
  end type

  ! Blocks
  type Cgns_Block_Type
    character(len=80)                     :: name
    integer                               :: type
    integer                               :: mesh_info(3)
    integer                               :: n_sects
    type(Cgns_Section_Type), allocatable  :: section(:)
    integer                               :: n_bnd_conds
    type(Cgns_Bnd_Cond_Type), allocatable :: bnd_cond(:)
    integer                               :: n_coords
    character(len=80)                     :: coord_name(3)
    integer                               :: n_solutions
    type(Cgns_Solution_Type), allocatable :: solution(:)
  end type

  ! Base
  integer :: n_bases
  type Cgns_Base_Type
    character(len=80)                  :: name
    integer                            :: cell_dim
    integer                            :: phys_dim
    integer                            :: n_blocks
    type(Cgns_Block_Type), allocatable :: block(:)
  end type
  type(Cgns_Base_Type), allocatable :: cgns_base(:)

  ! Some global counters (this is a bit ugly)
  integer :: cnt_nodes
  integer :: cnt_cells
  integer :: cnt_blocks     ! probably not needed
  integer :: cnt_bnd_cells

  integer :: cnt_hex
  integer :: cnt_pyr
  integer :: cnt_wed
  integer :: cnt_tet
  integer :: cnt_qua
  integer :: cnt_tri

  integer :: cnt_x
  integer :: cnt_y
  integer :: cnt_z

  integer, allocatable :: tflows_2_cgns_cells(:)

  ! Block-wise counter of boundary cells
  integer :: cnt_block_bnd_cells  ! probably not needed
  integer :: cnt_bnd_conds
  character(len=80) :: bnd_cond_names(1024)

  contains

  ! Common
  include 'Cgns_Mod/Initialize_Counters.f90'

  ! Seq only
  include 'Cgns_Mod/Read_Base_Info.f90'
  include 'Cgns_Mod/Read_Number_Of_Bases_In_File.f90'
  include 'Cgns_Mod/Read_Number_Of_Blocks_In_Base.f90'
  include 'Cgns_Mod/Read_Block_Info.f90'
  include 'Cgns_Mod/Read_Block_Type.f90'
  include 'Cgns_Mod/Read_Number_Of_Element_Sections.f90'
  include 'Cgns_Mod/Read_Section_Info.f90'
  include 'Cgns_Mod/Read_Number_Of_Bnd_Conds_In_Block.f90'
  include 'Cgns_Mod/Read_Bnd_Conds_Info.f90'
  include 'Cgns_Mod/Read_Number_Of_Coordinates_In_Block.f90'
  include 'Cgns_Mod/Read_Coordinate_Info.f90'
  include 'Cgns_Mod/Read_Coordinate_Array.f90'
  include 'Cgns_Mod/Read_Section_Connections.f90'
  include 'Cgns_Mod/Merge_Nodes.f90'

  ! Par only
  include 'Cgns_Mod/Get_Arrays_Dimensions_Par.f90'

  ! Par & Seq pairs

  include 'Cgns_Mod/Open_File_Seq.f90'
  include 'Cgns_Mod/Open_File_Par.f90'

  include 'Cgns_Mod/Close_File_Seq.f90'
  include 'Cgns_Mod/Close_File_Par.f90'

  include 'Cgns_Mod/Write_Base_Info_Seq.f90'
  include 'Cgns_Mod/Write_Base_Info_Par.f90'

  include 'Cgns_Mod/Write_Block_Info_Seq.f90'
  include 'Cgns_Mod/Write_Block_Info_Par.f90'

  include 'Cgns_Mod/Write_Coordinate_Array_Seq.f90'
  include 'Cgns_Mod/Write_Coordinate_Array_Par.f90'

  include 'Cgns_Mod/Write_Section_Connections_Seq.f90'
  include 'Cgns_Mod/Write_Section_Connections_Par.f90'

  include 'Cgns_Mod/Write_Solution_Info_Seq.f90'
  include 'Cgns_Mod/Write_Solution_Info_Par.f90'
  
  include 'Cgns_Mod/Write_Field_Seq.f90'
  include 'Cgns_Mod/Write_Field_Par.f90'


  ! move this two from Processor to Shared?
  !include 'Cgns_Mod/Save_Grid_Seq.f90'
  !include 'Cgns_Mod/Save_Grid_Par.f90'
  !include 'Cgns_Mod/Add_Fields_To_Grid_Seq.f90'
  !include 'Cgns_Mod/Add_Fields_To_Grid_Par.f90'
  end module
