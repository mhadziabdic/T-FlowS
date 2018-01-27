module cgns_mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Modules]-----------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!------------------------------------------------------------------------------!

  ! private
  !integer, private  :: 

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
  integer           :: first_node
  integer           :: last_node
  real, allocatable :: x_coord(:)
  real, allocatable :: y_coord(:)
  real, allocatable :: z_coord(:)

  ! elements
  integer              :: sect_id
  integer              :: n_sects
  integer              :: n_cells
  integer              :: first_cell
  integer              :: last_cell
  integer              :: cell_type
  integer              :: n_bnd
  character(len=11)    :: sect_name
  integer              :: iparent_flag
  integer              :: iparent_data
  integer, allocatable :: hexa_connections(:, :)
  integer, allocatable :: pyra_connections(:, :)
  integer, allocatable :: pris_connections(:, :)
  integer, allocatable :: tetr_connections(:, :)
  integer, allocatable :: tria_connections(:, :)
  integer, allocatable :: para_connections(:, :)

  ! bc
  integer              :: bc_idx
  integer              :: bc_type
  integer              :: ndataset
  integer              :: NormalListFlag
  integer              :: bt_set_type
  integer              :: n_nodes_in_bc
  integer              :: NormalIndex(3)
  integer, allocatable :: bc_points(:)
  integer              :: n_b_nodes_meth_1
  integer              :: n_b_nodes_meth_2
  integer              :: data_type
  integer              :: n_bc
  integer              :: bc_loc
  integer, allocatable :: bc_mark(:)
  character(len=11)    :: bc_name

  ! buffers
  real*4,  allocatable :: buffer_single(:)
  real*8,  allocatable :: buffer_double(:)
  integer, allocatable :: buffer_r2(:,:)
  integer, allocatable :: buffer_r1(:)
  
  integer              :: ier
  integer              :: i
  integer              :: j

  ! private

  ! functions
  public :: Cgns_Mod_Open_File
  public :: Cgns_Mod_Read_Number_Of_Bases_In_File
  public :: Cgns_Mod_Print_Base_Info
  public :: Cgns_Mod_Read_Number_Of_Zones_In_Base
  public :: Cgns_Mod_Read_Zone_Info
  public :: Cgns_Mod_Read_Zone_Type
  public :: Cgns_Mod_Read_Number_Of_Coordinates
  public :: Cgns_Mod_Print_Coordinate_Info
  public :: Cgns_Mod_Read_Coordinate_Array
  public :: Cgns_Mod_Read_Number_Of_Element_Sections

  public :: Cgns_Mod_Read_Section_Info
  public :: Cgns_Mod_Read_Section_Connections

  public ::Cgns_Mod_Read_BC_Number_In_Zone
  public ::Cgns_Mod_Read_BC_Info
  public ::Cgns_Mod_Read_BC

  !public :: 
  !public :: 

  public :: Cgns_Mod_Append_Real_Buf_To_Ar_R1
  public :: Cgns_Mod_Append_Int_Buf_To_Ar_R1
  public :: Cgns_Mod_Append_Int_Buf_To_Ar_R2
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

  print "(A,A)", "# Reading the file:", trim(file_name)

  ! Open a CGNS file 
  call Cg_Open_F(file_name,    &  ! file name
                 CG_MODE_READ, &  ! open for read
                 file_id,      &  ! cgns file index number
                 ier)             ! error status

  if (ier .ne. 0) then
    print "(A,A)", "# FAILED to read the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

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
      print "(A)", "# Failed to get bases number"
      call Cg_Error_Exit_F()
    endif

  print "(A,I5)", "# Number of bases: ", n_bases

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
      print "(A)", "#   FAILED to get base info"
      call Cg_Error_Exit_F()
    endif

  print "(A,A)",   "#   Base Name: ",      base_name
  print "(A,I1)",  "#   Cell dimension: ", cell_dim
  print "(A,I1)",  "#   Phys dimension: ", phys_dim 

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
    print "(A)", "#   FAILED to get zones number"
    call Cg_Error_Exit_F()
  endif

  print "(A,I3)","#   Total zones:", n_zones

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
    print "(A)", "#     Failed read zone info"
    call Cg_Error_Exit_F()
  endif

  n_nodes = mesh_info(1)
  n_cells = mesh_info(2)

  print "(A,I3)",  "#     Zone index: ",                zone_id
  print "(A,A)",   "#     Zone name: ",                 zone_name
  print "(A,I16)", "#     Nodes: ",                     n_nodes
  print "(A,I16)", "#     Cells: ",                     n_cells
  print "(A,I16)", "#     Boundary nodes(if sorted): ", mesh_info(3)

  if (mesh_info(3) .ne. 0) then
    print "(A)", "#     B.C. nodes != 0 -> Unsupported"
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
    print "(A)", "#     Failed to get zone type"
    call Cg_Error_Exit_F()
  endif

  print "(A,A)", "#     Zone type is ", ZoneTypeName(zone_type)

  if (zone_type .eq. Structured) then
    print "(A)", "#     Structured cgns meshed are unsupported"
    stop
  endif

  print "(A,I9)", "#     Index dimension = ", 1

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Number_Of_Coordinates
!------------------------------------------------------------------------------!
!   Gets n_zones from base node base_id
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
    print "(A)", "#       FAILED to get number of coordinate arrays"
    call Cg_Error_Exit_F()
  endif

  print "(A,I9)","#       Number of coordinate arrays for zone:", n_coords

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Print_Coordinate_Info
!------------------------------------------------------------------------------!
!   Gets n_zones from base node base_id
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  !Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,          & ! cgns file index number
                       base_id,          & ! base index number
                       zone_id,          & ! zone index number
                       coord_id,         & ! Coordinate array index number
                       coord_data_type,  & ! realsingle or realdouble
                       coord_name,       & ! name of the coordinate array
                       ier)                ! error status

  if (ier .ne. 0) then
    print "(A)", "#       FAILED to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  print "(A,I3)", "#       Coord. id:",         coord_id
  print "(A,A)",  "#       Coord. Data Type: ", DataTypeName(coord_data_type)
  print "(A,A)",  "#       Name: ",             coord_name

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Coordinate_Array
!------------------------------------------------------------------------------!
!   Read grid coordinates (RealSingle)                                         !
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  if (coord_data_type .eq. RealSingle) then

    ! Allocate buffer_single array
    if(.not. allocated(buffer_single)) then
      allocate(buffer_single(first_node:last_node))
    end if

    ! Read grid coordinates (RealSingle)
    call Cg_Coord_Read_F(file_id,       & ! cgns file index number
                         base_id,       & ! base index number
                         zone_id,       & ! zone index number
                         coord_name,    & ! name of the coordinate array
                         RealSingle,    & ! realsingle or realdouble
                         first_node,    & ! lower range index
                         last_node,     & ! upper range index
                         buffer_single, & ! array of coordinate values
                         ier)             ! error status

  else if (coord_data_type .eq. RealDouble) then

    ! Allocate buffer_double array
    if(.not. allocated(buffer_double)) &
      allocate(buffer_double(first_node:last_node))

    ! Read grid coordinates (RealDouble)
    call Cg_Coord_Read_F(file_id,       & ! cgns file index number
                         base_id,       & ! base index number
                         zone_id,       & ! zone index number
                         coord_name,    & ! name of the coordinate array
                         RealDouble,    & ! realsingle or realdouble
                         first_node,    & ! lower range index
                         last_node,     & ! upper range index
                         buffer_double, & ! array of coordinate values
                         ier)             ! error status

  end if

  if (ier.ne.0) then
    print *,"#       FAILED to read DoubleReal Coord"
    call Cg_Error_Exit_F()
  endif

  ! Append this data to coordinates arrays through buffer
  if (coord_data_type .eq. RealSingle) then
    allocate(buffer_double(first_node:last_node))

    buffer_double(first_node:last_node) &
      = real(buffer_single(first_node:last_node), 8)

    deallocate(buffer_single)
  end if

  select case ( coord_id )
   case (1)
      call Cgns_Mod_Append_Real_Buf_To_Ar_R1 (x_coord, buffer_double)
    case (2)
      call Cgns_Mod_Append_Real_Buf_To_Ar_R1 (y_coord, buffer_double)
    case (3)
      call Cgns_Mod_Append_Real_Buf_To_Ar_R1 (z_coord, buffer_double)
  end select

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

  print "(A)", "# ---------------------------"
  print "(A,I9)", "#       Number of sections: ", n_sects

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

  print "(A,I9)", "#         section idx:  ", n_sects
  print "(A,A)",  "#         section name: ", sect_name
  print "(A,A)",  "#         section type: ", ElementTypeName(cell_type)
  print "(A,I9)", "#         first idx:",     first_cell
  print "(A,I9)", "#         last idx:",      last_cell

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_Section_Connections
!------------------------------------------------------------------------------!
!   Read elements connection for current sect_id
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  ! Allocate correct size buffer for reading
  if      (ElementTypeName(cell_type) .eq. 'HEXA_8')  then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:8, first_cell:last_cell))
  else if (ElementTypeName(cell_type) .eq. 'PYRA_5')  then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:5, first_cell:last_cell))
  else if (ElementTypeName(cell_type) .eq. 'PENTA_6') then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:6, first_cell:last_cell))
  else if (ElementTypeName(cell_type) .eq. 'TETRA_4') then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:4, first_cell:last_cell))
  else if (ElementTypeName(cell_type) .eq. 'QUAD_4')  then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:4, first_cell:last_cell))
      n_b_nodes_meth_1 = n_b_nodes_meth_1 + last_cell - first_cell + 1
  else if (ElementTypeName(cell_type) .eq. 'TRI_3')   then
    if(.not.allocated(buffer_r2)) &
      allocate(buffer_r2(1:3, first_cell:last_cell))
      n_b_nodes_meth_1 = n_b_nodes_meth_1 + last_cell - first_cell + 1
  end if

  ! Read element data (HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
  call Cg_Elements_Read_F(file_id,      & ! cgns file index number
                          base_id,      & ! base index number
                          zone_id,      & ! zone index number
                          sect_id,      & ! element section index
                          buffer_r2,    & ! element connectivity data
                          iparent_data, & ! for boundary or interface 
                          ier)            ! error status

  ! Append cell connections from buffer to hex array
  if      (ElementTypeName(cell_type) .eq. 'HEXA_8')  then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (hexa_connections ,buffer_r2)
  else if (ElementTypeName(cell_type) .eq. 'PYRA_5')  then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (pyra_connections ,buffer_r2)
  else if (ElementTypeName(cell_type) .eq. 'PENTA_6') then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (pris_connections ,buffer_r2)
  else if (ElementTypeName(cell_type) .eq. 'TETRA_4') then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (tetr_connections ,buffer_r2)
  else if (ElementTypeName(cell_type) .eq. 'QUAD_4')  then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (para_connections ,buffer_r2)
  else if (ElementTypeName(cell_type) .eq. 'TRI_3')   then
    call Cgns_Mod_Append_Int_Buf_To_Ar_R2 (tria_connections ,buffer_r2)
  end if

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_BC_Number_In_Zone
!------------------------------------------------------------------------------!
!   Read B.C. info for zone_id
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  ! Get info for an element section
  ! Recieves sect_name, first_cell: last_cell, n_bnd and cell_type
  ! Get number of boundary condition in zone
    call Cg_Nbocos_F(file_id, & ! cgns file index number
                     base_id, & ! base index number
                     zone_id, & ! zone index number
                     n_bc,    & ! number of boundary conditions in zone
                     ier)       ! error status
    if (ier.ne.0) then
      print *,"#     FAILED to obtain number of b.c."
      call Cg_Error_Exit_F()
    endif


  ! Read grid location (where b.c. is placed)
  call Cg_Gridlocation_Read_F(bc_loc, & ! Location in the grid
                              ier)      ! error status

  ! I do not know why, but order of rows below matters
  print "(A,A)",  "#     B.C. locationn: ", GridLocationName(bc_loc)

  print "(A,I3)", "#     B.C. types: ", n_bc
  print "(A,I3)", "#     Zone index: ", zone_id

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_BC_Info
!------------------------------------------------------------------------------!
!   Read boundary condition info for zone_id
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  ! Get boundary condition info 
  call Cg_Boco_Info_F(file_id,        & ! cgns file index number
                      base_id,        & ! base index number
                      zone_id,        & ! zone index number
                      bc_idx,         & ! boundary condition index number
                      bc_name,        & ! name of the boundary condition
                      bc_type,        & ! type of boundary condition
                      bt_set_type,    & ! the extent of the boundary condition
                      n_nodes_in_bc,  & ! number of points/cells to make BC
                      NormalIndex,    & ! Index vector normal direction of BC
                      NormalListFlag, & ! flag indicating if the normals are defined in NormalList 
                      data_type,      & ! data type used in the definition of the normals
                      ndataset,       & ! number of boundary condition datasets
                      ier)              ! error status
  if (ier.ne.0) then
    print *,"#     FAILED to obtain number of b.c."
    call Cg_Error_Exit_F()
  endif


  print "(A,I3)", "#       B.C. index: ",            bc_idx
  print "(A,A)",  "#       B.C. name: ",             bc_name
  !print "(A,A)",  "#       B.C. type: ",             BCTypeName(bc_type)
  !print "(A,I3)", "#       B.C. nodes: ",            n_nodes_in_bc
  print "(A,A)",  "#       B.C. Extent: ",           PointSetTypeName(bt_set_type)
  !print "(A,9I9)","#       B.C. normal: ",           NormalIndex(1:3)
  !print "(A,5I5)","#       B.C. normal list flag: ", NormalListFlag
  !print "(A,A)",  "#       B.C. normal data_type: ", DataTypeName(data_type)
  !print "(A,I3)", "#       B.C. datasets: ",         ndataset

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Read_BC
!------------------------------------------------------------------------------!
!   Read b.c. for bc_id
!------------------------------------------------------------------------------!
  implicit none  
!------------------------------------------------------------------------------!

  ! B.C. list can be given by two type of array
  ! 1) Based on "PointRange" 
  ! ("PointRange constitutes rectangular region")
  ! 2) Based on "PointList"
  ! ("In all other cases, PointList")
  if     (bt_set_type .eq. PointRange) then
    ! For a bt_set_type of PointRange, n_nodes_in_bc is always 2
    if (.not. allocated(buffer_r2)) allocate(buffer_r2(2,1))

    ! Read boundary condition data and normals 
    call Cg_Boco_Read_F(file_id,        & ! cgns file index number
                        base_id,        & ! base index number
                        zone_id,        & ! zone index number
                        bc_idx,         & ! boundary condition index number
                        buffer_r2,      & ! array of point or element indices defining the boundary condition region
                        NormalListFlag, & ! list of vectors normal to the boundary condition patch pointing into the interior of the zone
                        ier)              ! error status

    first_node = buffer_r2(1,1)
    last_node  = buffer_r2(2,1)

    n_b_nodes_meth_2 = n_b_nodes_meth_2 + last_node - first_node + 1

    print "(A,I6,A,I6)", "#       Point Range: ", first_node, " -", last_node

    ! Buffer for bc_points(to do)
 

    ! Buffer for bc_mark
    if (.not. allocated(buffer_r1)) & 
      allocate(buffer_r1(first_node:last_node))

    first_cell = lbound(x_coord, dim=1)
    last_cell  = ubound(x_coord, dim=1)

    ! Color nodes with b.c.
    do i = first_node,  last_node - first_node + 1
      do j = first_cell, last_cell
        if ( j .eq. i) buffer_r1(j) = bc_idx
      end do
    end do
  
  else if (bt_set_type .eq. PointList) then
    ! For bt_set_type of PointList, n_nodes_in_bc is the number of points/ elements in the list
    if (.not. allocated(buffer_r1)) &
      allocate(buffer_r1(1:n_nodes_in_bc))

    ! Read boundary condition data and normals 
    call Cg_Boco_Read_F(file_id,       & ! cgns file index number
                        base_id,       & ! base index number
                        zone_id,       & ! zone index number
                        bc_idx,         & ! boundary condition index number
                        buffer_r1,      & ! array of point or element indices defining the boundary condition region
                        NormalListFlag, & ! list of vectors normal to the boundary condition patch pointing into the interior of the zone
                        ier)              ! error status

    n_b_nodes_meth_2 = n_b_nodes_meth_2 + n_nodes_in_bc
    print "(A,I6,A,I6)", "#       Point List s/f: ", 1, " -", n_nodes_in_bc

    ! Color nodes with b.c. (to do)


  end if

  !call Cgns_Mod_Append_Int_Buf_To_Ar_R1(bc_points, buffer_r1)
  call Cgns_Mod_Append_Int_Buf_To_Ar_R1(bc_mark,   buffer_r1)

  print "(A)", "# ---------------------------"

end
!------------------------------------------------------------------------------!









!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Append_Real_Buf_To_Ar_R1(ar_r1, buffer_r1)
!------------------------------------------------------------------------------!
!   Appends one real*8 array to another                                        !
!------------------------------------------------------------------------------!
  implicit none  
!---------------------------------[Arguments]----------------------------------!
  real*8, allocatable  :: ar_r1(:)
  real*8, allocatable  :: buffer_r1(:)
  real*8, allocatable  :: tmp_buffer_real(:)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: x1, x2, y1, y2
!------------------------------------------------------------------------------!

  y1 = lbound(buffer_r1, dim=1) ! buffer_r1
  y2 = ubound(buffer_r1, dim=1) ! buffer_r1

  if (.not.allocated(ar_r1)) then
    allocate(ar_r1(y1:y2))

    ar_r1(y1:y2) = buffer_r1(y1:y2)
  else
    x1 = lbound(ar_r1, dim=1) ! ar_r1
    x2 = ubound(ar_r1, dim=1) ! ar_r1

    allocate(tmp_buffer_real(x1:x2))
    tmp_buffer_real(x1:x2) = ar_r1(x1:x2)

    deallocate(ar_r1)
    allocate(ar_r1(x1:x2 + y2 - y1 + 1))
    
    ar_r1(x1:x2) = tmp_buffer_real(x1:x2)
    deallocate(tmp_buffer_real)

    ar_r1(x2+1:x2 + y2 - y1 + 1) = buffer_r1(y1:y2)
  end if

  deallocate(buffer_r1)

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Append_Int_Buf_To_Ar_R1(ar_r1, buffer_r1)
!------------------------------------------------------------------------------!
!   Appends one int array to another                                        !
!------------------------------------------------------------------------------!
  implicit none  
!---------------------------------[Arguments]----------------------------------!
  integer, allocatable  :: ar_r1(:)
  integer, allocatable  :: buffer_r1(:)
  integer, allocatable  :: tmp_buffer_real(:)
!-----------------------------------[Locals]-----------------------------------!
  integer               :: x1, x2, y1, y2
!------------------------------------------------------------------------------!

  y1 = lbound(buffer_r1, dim=1) ! buffer_r1
  y2 = ubound(buffer_r1, dim=1) ! buffer_r1


  if (.not.allocated(ar_r1)) then
    allocate(ar_r1(y1:y2))

    ar_r1(y1:y2) = buffer_r1(y1:y2)
  else
    x1 = lbound(ar_r1, dim=1) ! ar_r1
    x2 = ubound(ar_r1, dim=1) ! ar_r1

    allocate(tmp_buffer_real(x1:x2))
    tmp_buffer_real(x1:x2) = ar_r1(x1:x2)

    deallocate(ar_r1)
    allocate(ar_r1(x1:x2 + y2 - y1 + 1))
    
    ar_r1(x1:x2) = tmp_buffer_real(x1:x2)
    deallocate(tmp_buffer_real)

    ar_r1(x2+1:x2 + y2 - y1 + 1) = buffer_r1(y1:y2)
  end if

  deallocate(buffer_r1)

end
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
subroutine Cgns_Mod_Append_Int_Buf_To_Ar_R2(ar_r2, buffer_ar_r2)
!------------------------------------------------------------------------------!
!   Appends one real*8 array to another                                        !
!------------------------------------------------------------------------------!
  implicit none  
!---------------------------------[Arguments]----------------------------------!
  integer, allocatable :: ar_r2(:,:)
  integer, allocatable :: buffer_ar_r2(:,:)
  integer, allocatable :: tmp_buffer_int(:,:)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: x1, x2, y1, y2, d1, d2
!------------------------------------------------------------------------------!

  y1 = lbound(buffer_ar_r2, dim=2) ! buffer_ar_r2
  y2 = ubound(buffer_ar_r2, dim=2) ! buffer_ar_r2
  d1 = lbound(buffer_ar_r2, dim=1) ! ar_r2
  d2 = lbound(buffer_ar_r2, dim=1) ! ar_r2

  if (.not.allocated(ar_r2)) then
    allocate(ar_r2(d1:d2, y1:y2))

    ar_r2(d1:d2, y1:y2) = buffer_ar_r2(d1:d2, y1:y2)
  else
    x1 = lbound(ar_r2, dim=2) ! ar_r2
    x2 = ubound(ar_r2, dim=2) ! ar_r2

    allocate(tmp_buffer_int(d1:d2, x1:x2))
    tmp_buffer_int(d1:d2, x1:x2) = ar_r2(d1:d2, x1:x2)

    deallocate(ar_r2)
    allocate(ar_r2(d1:d2, x1:x2 + y2 - y1 + 1))
    
    ar_r2(d1:d2, x1:x2) = tmp_buffer_int(d1:d2, x1:x2)
    deallocate(tmp_buffer_int)

    ar_r2(d1:d2, x2+1:x2 + y2 - y1 + 1) = buffer_ar_r2(d1:d2, y1:y2)
  end if

  deallocate(buffer_ar_r2)

end
!------------------------------------------------------------------------------!

end module