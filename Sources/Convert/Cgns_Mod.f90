!==============================================================================!
  module Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Modules]-----------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!==============================================================================!

  ! file
  integer           :: file_id
  character(len=80) :: file_name

  ! base
  integer           :: base_id
  integer           :: n_bases
  character(len=80) :: base_name

  ! zone
  integer           :: zone_id
  integer           :: n_zones
  character(len=80) :: zone_name
  integer           :: zone_type
  integer           :: mesh_info(3)

  ! coordinates
  integer           :: n_coords
  integer           :: n_nodes
  integer           :: coord_id
  integer           :: coord_data_type
  character(len=80) :: coord_name
  real, allocatable :: x_coord(:)
  integer           :: last_x
  real, allocatable :: y_coord(:)
  integer           :: last_y
  real, allocatable :: z_coord(:)
  integer           :: last_z

  ! elements
  integer              :: sect_id
  integer              :: n_sects
  integer              :: n_cells
  integer              :: first_cell
  integer              :: last_cell
  integer              :: cell_type
  integer              :: n_bnd
  character(len=80)    :: sect_name
  integer              :: iparent_flag
  integer              :: iparent_data
  integer              :: n_hexa
  integer              :: last_hexa
  integer, allocatable :: hexa_connections(:, :)
  integer              :: n_pyra
  integer              :: last_pyra
  integer, allocatable :: pyra_connections(:, :)
  integer              :: n_pris
  integer              :: last_pris
  integer, allocatable :: pris_connections(:, :)
  integer              :: n_tetr
  integer              :: last_tetr
  integer, allocatable :: tetr_connections(:, :)
  integer              :: n_tria
  integer              :: last_tria
  integer, allocatable :: tria_connections(:, :)
  integer              :: n_quad
  integer              :: last_quad
  integer, allocatable :: quad_connections(:, :)

  ! bc
  integer              :: bc_id(50)
  integer, allocatable :: bc_mark(:)
  character(len=80)    :: bc_name

  ! buffers
  real   , allocatable :: buffer_double(:)
  integer, allocatable :: buffer_r2(:,:)

  integer              :: ier, c, i, j, k

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
