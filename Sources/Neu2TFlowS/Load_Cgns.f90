!==============================================================================!
  subroutine Load_Cgns(grid) 
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!   https://cgns.github.io/CGNS_docs_current/midlevel/structural.html
!   |-> mesh_info          !
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use neu_mod 
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!-----------------------------------[Locals]-----------------------------------!
  character(len=80)       :: name_in
  integer                 :: i, j, n_blocks, n_bnd_sect, dum1, dum2
  integer,allocatable     :: temp(:)
  integer                 :: c, n, dir
  integer                 :: file_idx, base_idx, zone_idx
  integer                 :: sect_idx
  integer                 :: ier
  integer                 :: n_sect
  integer                 :: n_bnd
  integer                 :: cell_type, iparent_data, iparent_flag

  character(len=80)       :: zone_name, name_sect
  character*32            :: base_name
  character*32            :: coord_name
  integer                 :: n_zones
  integer                 :: index_dim
  integer                 :: cell_dim
  integer                 :: phys_dim
  integer                 :: n_bases
  integer                 :: zone_type
  integer                 :: n_nodes
  integer                 :: n_cells
  integer                 :: n_b_cells_meth_1
  integer                 :: n_b_cells_meth_2
  integer                 :: n_bc
  integer                 :: n_coords
  integer                 :: data_type
  integer                 :: coord_idx
  integer                 :: last_node
  integer                 :: first_node

  ! *8
  integer               :: first_cell
  integer               :: last_cell
  integer               :: mesh_info(3)

  real*4,  allocatable :: buffer_single(:)
  real*8,  allocatable :: buffer_double(:)
  integer, allocatable :: buffer_r2(:,:)
  real*8,  allocatable :: x_array(:), y_array(:), z_array(:)

  integer, allocatable :: hex_connections(:, :)
  integer, allocatable :: pyramid_connections(:, :)
  integer, allocatable :: prism_connections(:, :)
  integer, allocatable :: tetra_connections(:, :)
  integer, allocatable :: tri_2d_connections(:, :)
  integer, allocatable :: box_2d_connections(:, :)
  integer, allocatable :: cell_connections(:, :)
  
  character*32 bc_name
  integer               :: bc_type, bc_idx
  integer ndataset
  integer NormalListFlag
  integer ptset_type
  integer n_nodes_in_bc
  integer NormalIndex(3)
  integer grid_loc
  integer, allocatable :: bc_points(:)
  integer, allocatable :: buffer_r1(:)
!==============================================================================!
!---------------------------------[Interface]----------------------------------!
  interface
    subroutine Append_Int_Ar_To_Array_Rank_1(data_ar_r1, data_ar_to_append_r1)
      implicit none  
      integer, allocatable :: data_ar_r1(:)
      integer, allocatable :: data_ar_to_append_r1(:)
      integer, allocatable :: tmp_buffer_real(:)
      integer              :: x1, x2, y1, y2
    end subroutine Append_Int_Ar_To_Array_Rank_1

    subroutine Append_Real_Ar_To_Array_Rank_1(data_ar_r1, data_ar_to_append_r1)
      implicit none  
      real*8, allocatable :: data_ar_r1(:)
      real*8, allocatable :: data_ar_to_append_r1(:)
      real*8, allocatable :: tmp_buffer_real(:)
      integer             :: x1, x2, y1, y2
    end subroutine Append_Real_Ar_To_Array_Rank_1

    subroutine Append_Int_Ar_To_Array_Rank_2(data_ar_r2, data_ar_to_append_r2)
      implicit none  
      integer, allocatable :: data_ar_r2(:,:)
      integer, allocatable :: data_ar_to_append_r2(:,:)
      integer, allocatable :: tmp_buffer_int(:,:)
      integer              :: x1, x2, y1, y2, d1, d2
    end subroutine Append_Int_Ar_To_Array_Rank_2

  end interface
!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = ".cgns"

  !Open a CGNS file 
  print "(A,A)", "# Reading the file:", trim(name_in)
  call Cg_Open_F(name_in,      &  ! file name
                 CG_MODE_READ, &  ! open for read
                 file_idx,     &  ! cgns file index number
                 ier)             ! error status
    if (ier .ne. 0) then
      print "(A,A)", "# FAILED to read the file: ", trim(name_in)
      call Cg_Error_Exit_F()
    endif
  
  ! Get number of CGNS base nodes in file 
  call Cg_Nbases_F(file_idx, & ! cgns file index number
                   n_bases,  & ! number of bases present in the CGNS file
                   ier)        ! error status
    if (ier .ne. 0) then
      print "(A)", "# Failed to get bases number"
      call Cg_Error_Exit_F()
    endif
  print "(A,I5)", "# Number of bases: ", n_bases

  n_nodes = 0 
  n_cells = 0 
  n_b_cells_meth_1 = 0 

  ! Browse through all bases
  bases: do base_idx = 1, n_bases
    print "(A)", "# ---------------------------"

    ! Read CGNS base information 
    call Cg_Base_Read_F(file_idx,  & ! cgns file index number
                        base_idx,  & ! base index number
                        base_name, & ! name of the base
                        cell_dim,  & ! cell dimensions (3->volume cell)
                        phys_dim,  & ! number of coordinates to create vector
                        ier)         ! error status
      if (ier .ne. 0) then
        print "(A)", "#   FAILED to get base info"
        call Cg_Error_Exit_F()
      endif
    print "(A,A)",  "#   Base Name: ",      base_name
    print "(A,I1)",  "#   Cell dimension: ", cell_dim
    print "(A,I1)",  "#   Phys dimension: ", phys_dim 

    ! Get number of zones in base
    call Cg_Nzones_F(file_idx, & ! cgns file index number
                     base_idx, & ! base index number
                     n_zones,  & ! number of zones present in base
                     ier)        ! error status
    if (ier .ne. 0) then
      print "(A)", "#   FAILED to get zones number"
      call Cg_Error_Exit_F()
    endif
    print "(A,I3)","#   Total zones:", n_zones

    ! Browse through all zones
    zones: do zone_idx = 1, n_zones
      print "(A)", "# ---------------------------"

      ! Read zone information
      call Cg_Zone_Read_F(file_idx,  & ! cgns file index number
                          base_idx,  & ! base index number
                          zone_idx,  & ! zone index number
                          zone_name, & ! name of the zone
                          mesh_info, & ! n_noded, n_cells, n_b_nodes(if sorted)
                          ier)         ! error status
        if (ier .ne. 0) then
          print "(A)", "#     Failed read zone info"
          call Cg_Error_Exit_F()
        endif
      print "(A,I3)",  "#     Zone index: ",                zone_idx
      print "(A,A)",   "#     Zone name: ",                 zone_name
      print "(A,I16)", "#     Nodes: ",                     mesh_info(1)
      print "(A,I16)", "#     Cells: ",                     mesh_info(2)
      print "(A,I16)", "#     Boundary nodes(if sorted): ", mesh_info(3)

      ! Total nodes and cells
      n_nodes    = n_nodes + mesh_info(1)
      n_cells    = n_cells + mesh_info(2)
      ! First and last nodes in this zone
      first_node = 1
      last_node  = mesh_info(1)

      ! Get type of zone (structured or unstructured) 
      call Cg_Zone_Type_F(file_idx,  & ! cgns file index number
                          base_idx,  & ! base index number
                          zone_idx,  & ! zone index number
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
      index_dim = 1

      !---------------------------!
      !   Read coordinates block  !
      !---------------------------!

      print "(A,I9)", "#     Index dimension = ", index_dim

      ! Get number of coordinate arrays
      call Cg_Ncoords_F(file_idx, & ! cgns file index number
                        base_idx, & ! base index number
                        zone_idx, & ! zone index number
                        n_coords, & ! number of coordinate arrays for zone
                        ier)        ! error status
        if (ier .ne. 0) then
          print "(A)", "#       FAILED to get number of coordinate arrays"
          call Cg_Error_Exit_F()
        endif
      print "(A,I9)","#       Number of coordinate arrays for zone:", n_coords

      ! Read x, y and z coordinates
      coordinates: do coord_idx = 1, n_coords
        print "(A)", "# ---------------------------"

        ! Get info about a coordinate array 
        call Cg_Coord_Info_F(file_idx,   & ! cgns file index number
                             base_idx,   & ! base index number
                             zone_idx,   & ! zone index number
                             coord_idx,  & ! Coordinate array index number
                             data_type,  & ! realsingle or realdouble
                             coord_name, & ! name of the coordinate array
                             ier)          ! error status

        print "(A,I3)", "#       Coord. idx:", coord_idx
        print "(A,A)",  "#       Data_type: ", DataTypeName(data_type)
        print "(A,A)",  "#       Name: ",      coord_name

        if (data_type .eq. RealSingle) then

          ! Allocate buffer_single array
          if(.not. allocated(buffer_single)) &
            allocate(buffer_single(first_node:last_node))

          ! Read grid coordinates (RealSingle)
          call Cg_Coord_Read_F(file_idx,      & ! cgns file index number
                               base_idx,      & ! base index number
                               zone_idx,      & ! zone index number
                               coord_name,    & ! name of the coordinate array
                               RealSingle,    & ! realsingle or realdouble
                               first_node,    & ! lower range index
                               last_node,     & ! upper range index
                               buffer_single, & ! array of coordinate values
                               ier)             ! error status
          if (ier.ne.0) then
            print *,"#       FAILED to read SingleReal Coord"
            call Cg_Error_Exit_F()
          endif

        elseif (data_type .eq. RealDouble) then
          ! Allocate buffer_double array
          if(.not. allocated(buffer_double)) &
            allocate(buffer_double(first_node:last_node))

          ! Read grid coordinates (RealDouble)
          call Cg_Coord_Read_F(file_idx,      & ! cgns file index number
                               base_idx,      & ! base index number
                               zone_idx,      & ! zone index number
                               coord_name,    & ! name of the coordinate array
                               RealDouble,    & ! realsingle or realdouble
                               first_node,    & ! lower range index
                               last_node,     & ! upper range index
                               buffer_double, & ! array of coordinate values
                               ier)             ! error status
          if (ier.ne.0) then
            print *,"#       FAILED to read DoubleReal Coord"
            call Cg_Error_Exit_F()
          endif

        end if

        ! Append this data to coordinates arrays through buffer
        if (data_type.eq.RealSingle) then
          allocate(buffer_double(first_node:last_node))
          buffer_double(:) = real(buffer_single(:), 8)
          deallocate(buffer_single)
        end if

        if ( coord_idx == 1) call Append_Real_Ar_To_Array_Rank_1 &
          (data_ar_r1 = x_array, data_ar_to_append_r1 = buffer_double)
        if ( coord_idx == 2) call Append_Real_Ar_To_Array_Rank_1 &
          (data_ar_r1 = y_array, data_ar_to_append_r1 = buffer_double)
        if ( coord_idx == 3) call Append_Real_Ar_To_Array_Rank_1 &
          (data_ar_r1 = z_array, data_ar_to_append_r1 = buffer_double)

      end do coordinates

      !---------------------!
      !   Read cells block  !
      !---------------------!

      ! Get number of element sections 
      call Cg_Nsections_F(file_idx, & ! cgns file index number
                          base_idx, & ! base index number
                          zone_idx, & ! zone index number
                          n_sect,   & ! number of element sections
                          ier)        ! error status
      print "(A)", "# ---------------------------"
      print "(A,I9)", "#       Number of sections: ", n_sect

      ! Browse through all sections to read elements
      elements_connection: do sect_idx = 1, n_sect

        ! Get info for an element section
        ! Recieves name_sect, first_cell: last_cell, n_bnd and cell_type
        call Cg_Section_Read_F(file_idx,     & ! cgns file index number
                               base_idx,     & ! base index number
                               zone_idx,     & ! zone index number
                               sect_idx,     & ! element section index
                               name_sect,    & ! name of the Elements_t node
                               cell_type,    & ! type of element
                               first_cell,   & ! index of first element
                               last_cell,    & ! index of last element
                               n_bnd,        & ! index of last boundary element
                               iparent_flag, & ! if the parent data are defined
                               ier)            ! error status

        print "(A,I9)", "#         section idx:  ", n_sect
        print "(A,A)",  "#         section name: ", name_sect
        print "(A,A)",  "#         section type: ", ElementTypeName(cell_type)
        print "(A,I9)", "#         first idx:",     first_cell
        print "(A,I9)", "#         last idx:",      last_cell

        ! Allocate correct size buffer for reading
        if     (ElementTypeName(cell_type) .eq. 'HEXA_8')  then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:8, first_cell:last_cell))
        elseif (ElementTypeName(cell_type) .eq. 'PYRA_5')  then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:5, first_cell:last_cell))
        elseif (ElementTypeName(cell_type) .eq. 'PENTA_6') then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:6, first_cell:last_cell))
        elseif (ElementTypeName(cell_type) .eq. 'TETRA_4') then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:4, first_cell:last_cell))
        elseif (ElementTypeName(cell_type) .eq. 'QUAD_4')  then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:4, first_cell:last_cell))
            n_b_cells_meth_1 = n_b_cells_meth_1 + last_cell - first_cell + 1
        elseif (ElementTypeName(cell_type) .eq. 'TRI_3')   then
          if(.not.allocated(buffer_r2)) &
            allocate(buffer_r2(1:3, first_cell:last_cell))
            n_b_cells_meth_1 = n_b_cells_meth_1 + last_cell - first_cell + 1
        end if

        ! Read element data (HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
        call Cg_Elements_Read_F(file_idx,     & ! cgns file index number
                                base_idx,     & ! base index number
                                zone_idx,     & ! zone index number
                                sect_idx,     & ! element section index
                                buffer_r2,   & ! element connectivity data
                                iparent_data, & ! for boundary or interface 
                                ier)            ! error status

        ! Append cell connections from buffer to hex array
        if     (ElementTypeName(cell_type) .eq. 'HEXA_8')  then
          call Append_Int_Ar_To_Array_Rank_2(hex_connections,buffer_r2)
        elseif (ElementTypeName(cell_type) .eq. 'PYRA_5')  then
          call Append_Int_Ar_To_Array_Rank_2(pyramid_connections,buffer_r2)
        elseif (ElementTypeName(cell_type) .eq. 'PENTA_6') then
          call Append_Int_Ar_To_Array_Rank_2(prism_connections,buffer_r2)
        elseif (ElementTypeName(cell_type) .eq. 'TETRA_4') then
          call Append_Int_Ar_To_Array_Rank_2(tetra_connections,buffer_r2)
        elseif (ElementTypeName(cell_type) .eq. 'QUAD_4')  then
          call Append_Int_Ar_To_Array_Rank_2(box_2d_connections,buffer_r2)
        elseif (ElementTypeName(cell_type) .eq. 'TRI_3')   then
          call Append_Int_Ar_To_Array_Rank_2(tri_2d_connections,buffer_r2)
        end if

        ! something is wrong here with yf17.cgns example:
        ! during sect_idx = 1, n_sect it changes sect_idx (debug this)
        ! Number of sections:        15
        ! section idx:         15
        ! section name: intake                                                                          
        ! section type: TRI_3                           
        ! first idx:   528916
        ! last idx:   528942
        ! Number of sections:     10372
        ! section idx:      10372
        ! section name: exhaust                                                                         
        ! section type: TRI_3                           
        ! first idx:   528943
        ! last idx:   528978
        print "(A)", "# ---------------------------"

      end do elements_connection
    
      !---------------!
      !   B.C. block  !
      !---------------!
      ! For pointwise this method duplicates info on b.c. obtained above
      ! But unfortunately for gridgen this is the only option
      ! to know b.c. nodes

      ! Access a node via label/name, index pairs 
      call Cg_Goto_F(file_idx,   & ! cgns file index number
                     base_idx,   & ! base index number
                     ier,        & ! error status
                     'Zone_t',   & ! node of zone type
                     zone_idx,   & ! zone index number
                     'ZoneBC_t', & ! search for node "ZoneBC_t"
                     1,          & ! ???
                     'end')        ! indicates end of call
      if (ier.ne.0) then
        print *,"#     FAILED to navigate ro ZoneBC node"
        call Cg_Error_Exit_F()
      endif
      
      ! Get number of boundary condition in zone
      call Cg_Nbocos_F(file_idx,  & ! cgns file index number
                       base_idx,  & ! base index number
                       zone_idx,  & ! zone index number
                       n_bc,      & ! number of boundary conditions in zone
                       ier)         ! error status
      if (ier.ne.0) then
        print *,"#     FAILED to obtain number of b.c."
        call Cg_Error_Exit_F()
      endif

      ! Read grid location 
      call Cg_Gridlocation_Read_F(grid_loc, & ! Location in the grid
                                  ier)        ! error status
        if (grid_loc .eq. FaceCenter) then
          print "(A)", "#     GridLocation refers to elements, not nodes"
        end if

      print "(A,I3)", "#     Zone index: ", zone_idx
      print "(A,I3)", "#     B.C. types: ", n_bc
      n_b_cells_meth_2 = 0

      boundary_conditions: do bc_idx = 1, n_bc

        ! Get boundary condition info 
        call Cg_Boco_Info_F(file_idx,       & ! cgns file index number
                            base_idx,       & ! base index number
                            zone_idx,       & ! zone index number
                            bc_idx,         & ! boundary condition index number
                            bc_name,        & ! name of the boundary condition
                            bc_type,        & ! type of boundary condition
                            ptset_type,     & ! the extent of the boundary condition
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

        n_b_cells_meth_2 = n_b_cells_meth_2 + n_nodes_in_bc

        print "(A,I3)", "#       B.C. index: ",            bc_idx
        print "(A,A)",  "#       B.C. name: ",             bc_name
        !print "(A,A)",  "#       B.C. type: ",             BCTypeName(bc_type)
        print "(A,I3)", "#       B.C. nodes: ",            n_nodes_in_bc
        !print "(A,A)",  "#       B.C. Extent: ",           PointSetTypeName(ptset_type)
        !print "(A,9I9)","#       B.C. normal: ",           NormalIndex(1:3)
        !print "(A,5I5)","#       B.C. normal list flag: ", NormalListFlag
        !print "(A,A)",  "#       B.C. normal data_type: ", DataTypeName(data_type)
        !print "(A,I3)", "#       B.C. datasets: ",         ndataset

        if (.not. allocated(buffer_r1)) allocate(buffer_r1(1:n_nodes_in_bc))
        call Append_Int_Ar_To_Array_Rank_1(bc_points,buffer_r1)
        
        ! Read boundary condition data and normals 
        call Cg_Boco_Read_F(file_idx,       & ! cgns file index number
                            base_idx,       & ! base index number
                            zone_idx,       & ! zone index number
                            bc_idx,         & ! boundary condition index number
                            buffer_r1,      & ! array of point or element indices defining the boundary condition region
                            NormalListFlag, & ! list of vectors normal to the boundary condition patch pointing into the interior of the zone
                            ier)              ! error status
        print *, "sdfsdf"

        print "(A)", "# ---------------------------"
      end do boundary_conditions

    end do zones

  end do bases

  print "(A,I9)", "#       Total    nodes read: ",  n_nodes 
  print "(A,I9)", "#       Total    cells read: ",  n_cells
  print "(A,I9)", "#       Total b. cells read: ",  n_b_cells_meth_1
  print "(A,I9)", "#       Total b. nodes read: ",  n_b_cells_meth_2
  grid % n_nodes     = n_nodes
  grid % n_cells     = n_cells
  grid % n_bnd_cells = n_b_cells_meth_1

  ! Allocate memory =--> carefull, there is no checking!
  ! Allocate memory =--> carefull, there is no checking!
  !call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  ! call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)


  !do c = 1, grid % n_cells
  !  print "(A)", "# Cell: ", c
  !  print "(A)", "# Nodes: ", (grid % cells_n(n,c), n=1,8)
  !end do

stop

  end subroutine
!cg_nsections - Get number of element sections
!cg_section_read - Get info for an element section
!cg_ElementDataSize - Get size of element connectivity data array
!cg_ElementPartialSize - Get size of element connectivity data array for partial read
!cg_elements_read - Read element data
!cg_elements_partial_read - Read subset of element data
!cg_npe - Get number of nodes for an element type 