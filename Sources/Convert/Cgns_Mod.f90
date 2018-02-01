module cgns_mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Modules]-----------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!------------------------------------------------------------------------------!

  ! file
  integer           :: file_id
  character(len=50) :: file_name

  ! base
  integer           :: base_id
  integer           :: n_bases
  character(len=50) :: base_name

  ! zone
  integer           :: zone_id
  integer           :: n_zones
  character(len=50) :: zone_name
  integer           :: zone_type
  integer           :: mesh_info(3)

  ! coordinates
  integer           :: n_coords
  integer           :: n_nodes
  integer           :: coord_id
  integer           :: coord_data_type
  character(len=11) :: coord_name
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
  character(len=40)    :: sect_name
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
  character(len=11)    :: bc_name

  ! buffers
  real   , allocatable :: buffer_double(:)
  integer, allocatable :: buffer_r2(:,:)

  integer              :: ier, c, i, j, k

!------------------------------------------------------------------------------!
contains
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Open_File
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  print *, "# Reading the file:", trim(file_name)

  ! Open a CGNS file
  call Cg_Open_F(file_name,    &  ! file name
                 CG_MODE_READ, &  ! open for read
                 file_id,      &  ! cgns file index number
                 ier)             ! error status

  if (ier .ne. 0) then
    print *, "# FAILED to read the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

  !  set initial values
  n_nodes = 0
  n_cells = 0
  n_hexa = 0
  n_pyra = 0
  n_pris = 0
  n_tetr = 0
  n_tria = 0
  n_quad = 0

  last_x = 0
  last_y = 0
  last_z = 0

  last_hexa = 0
  last_pyra = 0
  last_pris = 0
  last_tetr = 0
  last_tria = 0
  last_quad = 0

  bc_id = 0

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Number_Of_Bases_In_File
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Get number of CGNS bases in file
  call Cg_Nbases_F(file_id, & ! cgns file index number
                   n_bases, & ! number of bases present in the CGNS file
                   ier)       ! error status

    if (ier .ne. 0) then
      print *, "# Failed to get bases number"
      call Cg_Error_Exit_F()
    endif

  print *, "# Number of bases: ", n_bases

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Print_Base_Info
!------------------------------------------------------------------------------!
!   Reads main info from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=40) :: base_name
  integer           :: phys_dim
  integer           :: cell_dim
!------------------------------------------------------------------------------!

  ! Read CGNS base information
  call Cg_Base_Read_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      base_name, & ! name of the base
                      cell_dim,         & ! cell dimensions (3->volume cell)
                      phys_dim,         & ! number of coordinates to create vector
                      ier)                ! error status

    if (ier .ne. 0) then
      print *, "#   FAILED to get base info"
      call Cg_Error_Exit_F()
    endif

  print *, "#   Base Name: ",      base_name
  print *, "#   Cell dimension: ", cell_dim
  print *, "#   Phys dimension: ", phys_dim

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Number_Of_Zones_In_Base
!------------------------------------------------------------------------------!
!   Gets n_zones from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

 ! Get number of zones in base
  call Cg_Nzones_F(file_id, & ! cgns file index number
                   base_id, & ! base index number
                   n_zones, & ! number of zones present in base
                   ier)       ! error status
  if (ier .ne. 0) then
    print *, "#   FAILED to get zones number"
    call Cg_Error_Exit_F()
  endif

  print *, "#   Total zones:", n_zones

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Zone_Info
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Read zone information
  call Cg_Zone_Read_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      zone_id,   & ! zone index number
                      zone_name, & ! name of the zone
                      mesh_info, & ! n_nodes, n_cells, n_b_nodes(if sorted)
                      ier)         ! error status

  if (ier .ne. 0) then
    print *, "#     Failed read zone info"
    call Cg_Error_Exit_F()
  endif

  ! total nodes and cells
  n_nodes = n_nodes + mesh_info(1)
  n_cells = n_cells + mesh_info(2)

  print *, "#     Zone index: ",                zone_id
  print *, "#     Zone name: ",                 zone_name
  print *, "#     Nodes: ",                     mesh_info(1)
  print *, "#     Cells: ",                     mesh_info(2)
  print *, "#     Boundary nodes(if sorted): ", mesh_info(3)

  if (mesh_info(3) .ne. 0) then
    print *, "#     B.C. nodes != 0 -> Unsupported"
    stop
  endif

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Zone_Type
!------------------------------------------------------------------------------!
!   Prints zone type                                                             !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Read type of zone
  call Cg_Zone_Type_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      zone_id,   & ! zone index number
                      zone_type, & ! structured or unstructured
                      ier)         ! error status
  if (ier .ne. 0) then
    print *, "#     Failed to get zone type"
    call Cg_Error_Exit_F()
  endif

  print *, "#     Zone type is ", ZoneTypeName(zone_type)

  if (zone_type .eq. Structured) then
    print *, "#     Structured cgns meshed are unsupported"
    stop
  endif

  print *, "#     Index dimension = ", 1

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Number_Of_Coordinates
!------------------------------------------------------------------------------!
!   Reads number of coordinates arrays from zone_id
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Get number of coordinate arrays (1 for unstructure, 3 for structured)
  call Cg_Ncoords_F(file_id,  & ! cgns file index number
                    base_id,  & ! base index number
                    zone_id,  & ! zone index number
                    n_coords, & ! number of coordinate arrays for zone
                    ier)        ! error status

  if (ier .ne. 0) then
    print *, "#       FAILED to get number of coordinate arrays"
    call Cg_Error_Exit_F()
  endif

  print *, "#       Number of coordinate arrays for zone:", n_coords

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Print_Coordinate_Info
!------------------------------------------------------------------------------!
!   Reads coord_name of coord_id
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,         & ! cgns file index number
                       base_id,         & ! base index number
                       zone_id,         & ! zone index number
                       coord_id,        & ! Coordinate array index number
                       coord_data_type, & ! realsingle or realdouble
                       coord_name,      & ! name of the coordinate array
                       ier)               ! error status

  if (ier .ne. 0) then
    print *, "#       FAILED to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  print *, "#       Coord. id:",         coord_id
  print *, "#       Coord. Data Type: ", DataTypeName(coord_data_type)
  print *, "#       Name: ",             coord_name

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Coordinate_Array
!------------------------------------------------------------------------------!
!   Read grid coordinates (RealDouble)                                         !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  i = 1
  j = mesh_info(1)

  allocate(buffer_double(i:j))

  ! Read grid x coordinates
  call Cg_Coord_Read_F(file_id,       & ! cgns file index number
                       base_id,       & ! base index number
                       zone_id,       & ! zone index number
                       coord_name,    & ! name of the coordinate array
                       RealDouble,    & ! realsingle or realdouble
                       i,             & ! lower range index
                       j,             & ! upper range index
                       buffer_double, & ! array of coordinate values
                       ier)           ! error status

  if (ier.ne.0) then
    print *, "#       FAILED to read DoubleReal Coord", coord_name
    call Cg_Error_Exit_F()
  endif

  select case (coord_id)
    case (1)
      i = last_x + 1
      j = last_x + j
      x_coord(i:j) = buffer_double(:)
      last_x = j
    case (2)
      i = last_y + 1
      j = last_y + j
      y_coord(i:j) = buffer_double(:)
      last_y = j
    case (3)
      i = last_z + 1
      j = last_z + j
      z_coord(i:j) = buffer_double(:)
      last_z = j
  end select

  deallocate(buffer_double)
end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Number_Of_Element_Sections
!------------------------------------------------------------------------------!
!   Gets n_sects from zone
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Get number of element sections
  call Cg_Nsections_F(file_id,  & ! cgns file index number
                      base_id,  & ! base index number
                      zone_id,  & ! zone index number
                      n_sects,  & ! number of element sections
                      ier)        ! error status

  if (ier.ne.0) then
    print *, "#     FAILED to read number of elements"
    call Cg_Error_Exit_F()
  endif

  print *, "# ---------------------------"
  print *, "#       Number of sections: ", n_sects

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Section_Info
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect_id
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Get info for an element section
  ! Recieves sect_name, first_cell: last_cell, n_bnd and cell_type
  call Cg_Section_Read_F(file_id,      & ! cgns file index number
                         base_id,      & ! base index number
                         zone_id,      & ! zone index number
                         sect_id,      & ! element section index
                         sect_name,    & ! name of the Elements_t node
                         cell_type,    & ! type of element
                         first_cell,   & ! index of first element
                         last_cell,    & ! index of last element
                         n_bnd,        & ! index of last boundary element
                         iparent_flag, & ! if the parent data are defined
                         ier)            ! error status

  if (ier.ne.0) then
    print *, "#     FAILED to read section ", sect_id, " info"
    call Cg_Error_Exit_F()
  endif

  print *, "#         section idx:  ", sect_id
  print *, "#         section name: ", sect_name
  print *, "#         section type: ", ElementTypeName(cell_type)
  print *, "#         first cell:",    first_cell
  print *, "#         last cell:",     last_cell

  i = last_cell - first_cell + 1 ! cells in this sections

  ! count cells in sect_id
  if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_hexa = n_hexa + i
  if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_pyra = n_pyra + i
  if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_pris = n_pris + i
  if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_tetr = n_tetr + i
  if ( ElementTypeName(cell_type) .eq. 'QUAD_4' ) then
    n_quad = n_quad + i
    print *, "#         This section was identified as b.c."
    bc_id(sect_id) = 1
  end if
  if ( ElementTypeName(cell_type) .eq. 'TRI_3'  ) then
    n_tria = n_tria + i
    print *, "#         This section was identified as b.c."
    bc_id(sect_id) = 1
  end if

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Section_Connections
!------------------------------------------------------------------------------!
!   Read elements connection for current sect_id
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  ! Get info for an element section
  ! Recieves sect_name, first_cell: last_cell, n_bnd and cell_type
  call Cg_Section_Read_F(file_id,      & ! cgns file index number
                         base_id,      & ! base index number
                         zone_id,      & ! zone index number
                         sect_id,      & ! element section index
                         sect_name,    & ! name of the Elements_t node
                         cell_type,    & ! type of element
                         first_cell,   & ! index of first element
                         last_cell,    & ! index of last element
                         n_bnd,        & ! index of last boundary element
                         iparent_flag, & ! if the parent data are defined
                         ier)            ! error status

  if (ier.ne.0) then
    print *, "#     FAILED to read info for section ", sect_id
    call Cg_Error_Exit_F()
  endif

  i = first_cell
  j = last_cell

  !print *, "#     ", file_id, base_id, zone_id, sect_id, sect_name,first_cell, last_cell, n_bnd, iparent_flag, ier

  ! Read cells
  if (ElementTypeName(cell_type) .eq. 'HEXA_8')  allocate(buffer_r2(1:8, i:j))
  if (ElementTypeName(cell_type) .eq. 'PYRA_5')  allocate(buffer_r2(1:5, i:j))
  if (ElementTypeName(cell_type) .eq. 'PENTA_6') allocate(buffer_r2(1:6, i:j))
  if (ElementTypeName(cell_type) .eq. 'TETRA_4') allocate(buffer_r2(1:4, i:j))
  if (ElementTypeName(cell_type) .eq. 'QUAD_4')  allocate(buffer_r2(1:4, i:j))
  if (ElementTypeName(cell_type) .eq. 'TRI_3')   allocate(buffer_r2(1:3, i:j))

  ! Read HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3 elements
  call Cg_Elements_Read_F( file_id,      & ! cgns file index number
                           base_id,      & ! base index number
                           zone_id,      & ! zone index number
                           sect_id,      & ! element section index
                           buffer_r2,    & ! element connectivity data
                           iparent_data, & ! for boundary or interface
                           ier)            ! error status

  if (ier.ne.0) then
    print *, "#     FAILED to read elemets in section ", sect_id
    call Cg_Error_Exit_F()
  endif

  if      (ElementTypeName(cell_type) .eq. 'HEXA_8') then
    j = last_hexa + j - i + 1
    i = last_hexa + 1
    hexa_connections(1:8, i:j) = buffer_r2(1:8, first_cell:last_cell) + last_hexa
    last_hexa = j
  else if (ElementTypeName(cell_type) .eq. 'PYRA_5') then
    j = last_pyra + j - i + 1
    i = last_pyra + 1
    pyra_connections(1:5, i:j) = buffer_r2(1:5, first_cell:last_cell) + last_pyra
    last_pyra = j
  else if (ElementTypeName(cell_type) .eq. 'PENTA_6') then
    j = last_pris + j - i + 1
    i = last_pris + 1
    pris_connections(1:6, i:j) = buffer_r2(1:6, first_cell:last_cell) + last_pris
    last_pris = j
  else if (ElementTypeName(cell_type) .eq. 'TETRA_4') then
    j = last_tetr + j - i + 1
    i = last_tetr + 1
    tetr_connections(1:4, i:j) = buffer_r2(1:4, first_cell:last_cell) + last_tetr
    last_tetr = j
  else if (ElementTypeName(cell_type) .eq. 'QUAD_4') then
    j = last_quad + j - i + 1
    i = last_quad + 1
    quad_connections(1:4, i:j) = buffer_r2(1:4, first_cell:last_cell) + last_quad
    last_quad = j
  else if (ElementTypeName(cell_type) .eq. 'TRI_3') then
    j = last_tria + j - i + 1
    i = last_tria + 1
    tria_connections(1:3, i:j) = buffer_r2(1:3, first_cell:last_cell) + last_tria
    last_tria = j
  end if

  !do c = first_cell, last_cell
  !  print "(A,8I9)", "buffer=", buffer_r2(:,c)
  !end do

  print *, i, ":", j, " copyed from", first_cell, ":", last_cell," type = ",ElementTypeName(cell_type)



  deallocate(buffer_r2)

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Mark_Bound_Cond
!------------------------------------------------------------------------------!
!   Mark nodes in x,y,z arrays from tria_connections & quad_connections        !
!   with sect_id                                                                !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  bc_mark = 0

  ! Color quad nodes with b.c.
  if ( bc_id(sect_id) .ne. 0) then
    do c = 1, n_nodes
      do i = 1, last_quad
        do j = 1, 4
          if (c .eq. quad_connections(j,i)) then
            bc_mark(c) = sect_id
          end if
        end do
      end do
    end do
  end if

  ! Color tria nodes with b.c.
  if ( bc_id(sect_id) .ne. 0) then
    do c = 1, n_nodes
      do i = 1, last_tria
        do j = 1, 3
          if (c .eq. quad_connections(j,i)) then
            bc_mark(c) = sect_id
          end if
        end do
      end do
    end do
  end if

end
!------------------------------------------------------------------------------!
end module