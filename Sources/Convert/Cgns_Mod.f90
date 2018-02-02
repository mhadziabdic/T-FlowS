!==============================================================================!
  module Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!==============================================================================!

  ! file
  integer*8         :: file_id
  character(len=80) :: file_name

  ! base
  integer*8         :: base_id
  integer*8         :: n_bases
  character(len=80) :: base_name

  ! zone
  integer*8         :: zone_id
  integer*8         :: n_zones
  character(len=80) :: zone_name
  integer*8         :: zone_type
  integer*8         :: mesh_info(3)

  ! coordinates
  integer*8         :: n_coords
  integer*8         :: n_nodes
  integer*8         :: coord_id
  integer*8         :: coord_data_type
  character(len=80) :: coord_name
  real, allocatable :: x_coord(:)
  integer*8         :: last_x
  real, allocatable :: y_coord(:)
  integer*8         :: last_y
  real, allocatable :: z_coord(:)
  integer*8         :: last_z

  ! elements
  integer*8              :: sect_id
  integer*8              :: n_sects
  integer*8              :: n_cells
  integer*8              :: first_cell
  integer*8              :: last_cell
  integer*8              :: cell_type
  integer*8              :: n_bnd
  character(len=80)      :: sect_name
  integer*8              :: iparent_flag
  integer*8              :: iparent_data
  integer*8              :: n_hexa
  integer*8              :: last_hexa
  integer*8, allocatable :: hexa_connections(:, :)
  integer*8              :: n_pyra
  integer*8              :: last_pyra
  integer*8, allocatable :: pyra_connections(:, :)
  integer*8              :: n_pris
  integer*8              :: last_pris
  integer*8, allocatable :: pris_connections(:, :)
  integer*8              :: n_tetr
  integer*8              :: last_tetr
  integer*8, allocatable :: tetr_connections(:, :)
  integer*8              :: n_tria
  integer*8              :: last_tria
  integer*8, allocatable :: tria_connections(:, :)
  integer*8              :: n_quad
  integer*8              :: last_quad
  integer*8, allocatable :: quad_connections(:, :)

  ! bc
  integer*8              :: bc_id(50)
  integer*8, allocatable :: bc_mark(:)
  character(len=80)      :: bc_name

  ! buffers
  real             , allocatable :: buffer_double(:)
  integer*8        , allocatable :: buffer_r2(:,:)

  integer*8            :: ier, c, k
  integer*8            :: i, j

  contains

  include 'Cgns_Mod/Mark_Bound_Cond.f90'
  include 'Cgns_Mod/Open_File.f90'
  include 'Cgns_Mod/Print_Base_Info.f90'
  include 'Cgns_Mod/Print_Coordinate_Info.f90'
  include 'Cgns_Mod/Read_Coordinate_Array.f90'
  include 'Cgns_Mod/Read_Number_Of_Bases_In_File.f90'
  include 'Cgns_Mod/Read_Number_Of_Coordinates.f90'
  include 'Cgns_Mod/Read_Number_Of_Element_Sections.f90'
  include 'Cgns_Mod/Read_Number_Of_Zones_In_Base.f90'
  include 'Cgns_Mod/Read_Section_Connections.f90'
  include 'Cgns_Mod/Read_Section_Info.f90'
  include 'Cgns_Mod/Read_Zone_Info.f90'
  include 'Cgns_Mod/Read_Zone_Type.f90'

  end module
