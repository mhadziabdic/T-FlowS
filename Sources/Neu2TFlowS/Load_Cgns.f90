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
  use cgns_mod 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(len=80)       :: name_in
  integer                 :: n_blocks, n_bnd_sect, dum1, dum2
  integer,allocatable     :: temp(:)
  integer                 :: c, n

  integer                 :: index_dim
  integer                 :: n_nodes
  integer                 :: n_cells

  integer, allocatable :: cell_connections(:, :)

  integer dir
  integer n_faces
  integer fn(6,4)
  logical is_on_face
  
!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = ".cgns"


  file_name = name_in

  ! Open a CGNS file 
  call Cgns_Mod_Open_File
  
  ! Get number of CGNS bases in file 
  call Cgns_Mod_Read_Number_Of_Bases_In_File

  n_nodes = 0 
  n_cells = 0 
  n_b_nodes_meth_1 = 0 

  !------------------------------!
  !   Browse through all bases   !
  !------------------------------!
  bases: &
  do base_id = 1, n_bases

    ! Read CGNS base information 
    call Cgns_Mod_Read_Base_Info

    ! Get number of zones in base
    call Cgns_Mod_Read_Number_Of_Zones_In_Base

    !------------------------------!
    !   Browse through all zones   !
    !------------------------------!
    zones: &
    do zone_id = 1, n_zones

      call Cgns_Mod_Read_Zone_Info

      ! Total nodes and cells
      n_nodes = n_nodes + mesh_info(1)
      n_cells = n_cells + mesh_info(2)

      ! First and last nodes in this zone
      first_node = 1
      last_node  = mesh_info(1)

      ! Read type of zone (structured or unstructured) 
      call Cgns_Mod_Read_Zone_Type

      !---------------------------!
      !   Read coordinates block  !
      !---------------------------!

      ! Get number of coordinate arrays (1 for unstructured, 3 for structured)
      call Cgns_Mod_Read_Number_Of_Coordinates

      ! Read x, y and z coordinates
      coordinates: &
      do coord_id = 1, n_coords

        !Get info about coordinate coord_id
        call Cgns_Mod_Read_Coordinate_Info

        ! Read grid coordinates (fill: x_, y_, z_array)
        call Cgns_Mod_Read_Coordinate_Array

      end do coordinates

      !---------------------!
      !   Read cells block  !
      !---------------------!

      ! Get number of element sections 
      call Cgns_Mod_Read_Number_Of_Element_Sections

      ! Browse through all sections to read elements
      elements_connection: &
      do sect_id = 1, n_sects

        ! Get info for an element section
        call Cgns_Mod_Read_Section_Info

        ! Read element data (HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
        call Cgns_Mod_Read_Section_Connections

        ! something is wrong here with yf17.cgns example:
        ! during sect_id = 1, n_sect it changes sect_id (debug this)
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

      end do elements_connection
    
      !---------------!
      !   B.C. block  !
      !---------------!
      ! For pointwise this method duplicates info on b.c. obtained above
      ! But unfortunately for gridgen this is the only option
      ! to know b.c. nodes

      ! not needed ? 

      !! Access a node via label/name, index pairs 
      !call Cg_Goto_F(file_id,   & ! cgns file index number
      !               base_id,   & ! base index number
      !               ier,        & ! error status
      !               'Zone_t',   & ! node of zone type
      !               zone_id,   & ! zone index number
      !               'ZoneBC_t', & ! search for node "ZoneBC_t"
      !               1,          & ! ???
      !               'end')        ! indicates end of call
      !if (ier.ne.0) then
      !  print *,"#     FAILED to navigate ro ZoneBC node"
      !  call Cg_Error_Exit_F()
      !endif
      
      ! Get number of boundary condition in zone
      call Cgns_Mod_Read_BC_Number_In_Zone

      boundary_conditions: &
      do bc_idx = 1, n_bc

        ! Get boundary condition info 
        call Cgns_Mod_Read_BC_Info

        call Cgns_Mod_Read_BC

      end do boundary_conditions

        ! print b.c. nodes
        !do j = first_node, last_node
        !  if (bc_mark(j) .ne. 0) print *, "bc_mark=", bc_mark(j), &
        !  "x=", x_coord(j), "y=", y_coord(j), "z=", z_coord(j)
        !end do


        !--------------------------------------!
        !  Conversion to T-FlowS B.C. format   !
        !--------------------------------------!

        ! Reorder elements connectivity according to neu_mod

        ! HEXA_8 (-> f8n)
        do c = lbound(hex_connections,dim=2), &
          ubound(hex_connections,dim=2) - lbound(hex_connections,dim=2) + 1

          i = hex_connections(4, c)
          hex_connections(4, c) = hex_connections(3, c)
          hex_connections(3, c) = i
          i = hex_connections(8, c)
          hex_connections(8, c) = hex_connections(7, c)
          hex_connections(7, c) = i
        end do

        ! PYRA_5 (-> f5n)
        do c =  lbound(pyramid_connections,dim=2), &
          ubound(pyramid_connections,dim=2) - lbound(pyramid_connections,dim=2) + 1

          i = pyramid_connections(4, c)
          pyramid_connections(4, c) = pyramid_connections(3, c)
          pyramid_connections(3, c) = i
        end do

        ! PENTA_6 (-> f6n)
        ! no changes

        ! TETRA_4 (-> f4n)
        ! no changes



!        fn = 0
!
!        do i = 1, n_cells
!
!          n_faces = 6
!          if( fn(6,1) == -1 ) n_faces = 5
!          if( fn(5,1) == -1 ) n_faces = 4
!
!          do dir = 1, n_faces
!
!            n_nodes = 4
!            if( fn(dir, 4) == -1 ) n_nodes = 3
!
!            is_on_face = .true.
!            do n = 1, n_nodes 
!              if( bc_mark( hex_connections( fn(dir, n), n  )) .ne. bc_idx ) is_on_face = .false.
!            end do
!     
!            if( is_on_face ) then
!              BCtype(i, dir) = bc_idx
!            end if   
!          end do
!        end do


    end do zones

  end do bases

  print "(A,I9)", "#       Total    nodes read: ",            n_nodes 
  print "(A,I9)", "#       Total    cells read: ",            n_cells
  print "(A,I9)", "#       Total b. nodes read(method 1): ",  n_b_nodes_meth_1
  print "(A,I9)", "#       Total b. nodes read(method 2): ",  n_b_nodes_meth_2
  grid % n_nodes     = n_nodes
  grid % n_cells     = n_cells
  grid % n_bnd_cells = n_b_nodes_meth_2

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
!cg_ElementDataSize - Get size of element connectivity data array
!cg_ElementPartialSize - Get size of element connectivity data array for partial read
!cg_elements_read - Read element data
!cg_elements_partial_read - Read subset of element data
!cg_npe - Get number of nodes for an element type 