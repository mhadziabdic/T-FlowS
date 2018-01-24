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
  integer                 :: file_idx, base_idx, zone_idx, sect_idx, ier
  integer                 :: n_sect, n_bnd
  integer                 :: first_cell, last_cell, itype, iparent_data, iparent_flag

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
  integer                 :: n_bound
  integer                 :: n_coords
  integer                 :: data_type
  integer                 :: coord_idx
  integer                 :: last_node
  integer                 :: first_node
  integer                 :: mesh_info(3)
  real*4, allocatable     :: data_single(:)
  real*8, allocatable     :: data_double(:)
!==============================================================================!

  name_in = problem_name
  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = ".cgns"

  !Open a CGNS file 
  write(*,"(A,A)") "# Reading the file:", trim(name_in)
  call Cg_Open_F(name_in,      &  ! file name
                 CG_MODE_READ, &  ! open for read
                 file_idx,     &  ! CGNS file index number
                 ier)             ! Error status
    if (ier .ne. 0) then
      write(*,"(A,A)") "# FAILED to read the file: ", trim(name_in)
      call Cg_Error_Exit_f()
    endif

  ! Get number of CGNS base nodes in file 
  call Cg_Nbases_F(file_idx, & ! CGNS file index number
                   n_bases,  & ! Number of bases present in the CGNS file
                   ier)        ! Error status
    if (ier .ne. 0) then
      write(*,"(A)") "# Failed to get bases number"
      call Cg_Error_Exit_f()
    endif
  write(*,"(A,I5)") "# Number of bases: ", n_bases

  n_nodes = 0 
  n_cells = 0 
  n_bound = 0 

  ! Browse through all bases
  bases: do base_idx = 1, n_bases
    write(*,"(A)") "# ---------------------------"

    ! Read CGNS base information 
    call Cg_Base_Read_F(file_idx, & ! CGNS file index number
                        base_idx, & ! Base index number
                        base_name, & ! Name of the base
                        cell_dim, & ! Cell dimensions (3->volume cell)
                        phys_dim, & ! Number of coordinates to create vector
                        ier)        ! Error status
      if (ier .ne. 0) then
        write(*,"(A)") "#   FAILED to get base info"
        call Cg_Error_Exit_f()
      endif
    write(*,"(A,A)")  "#   Base Name: ",      base_name
    write(*,"(A,I1)") "#   Cell dimension: ", cell_dim
    write(*,"(A,I1)") "#   Phys dimension: ", phys_dim 

    ! Get number of zones in base
    call Cg_Nzones_F(file_idx, & ! CGNS file index number
                     base_idx, & ! Base index number
                     n_zones,  & ! Number of zones present in base
                     ier)        ! Error status
    if (ier .ne. 0) then
      write(*,"(A)") "#   FAILED to get zones number"
      call Cg_Error_Exit_f()
    endif
    write(*,"(A,I3)")"#   Total zones:", n_zones

    ! Browse through all zones
    zones: do zone_idx = 1, n_zones
      write(*,"(A)") "# ---------------------------"

      ! Read zone information
      call Cg_Zone_Read_F(file_idx,  & ! CGNS file index number
                          base_idx,  & ! Base index number
                          zone_idx,  & ! Zone index number
                          zone_name, & ! Name of the zone
                          mesh_info, & ! n_noded, n_cells, n_b_nodes(if sorted)
                          ier)         ! Error status
        if (ier .ne. 0) then
          write(*,"(A)") "#     Failed read zone info"
          call Cg_Error_Exit_f()
        endif
      write(*,"(A,I3)") "#     Zone index: ",                zone_idx
      write(*,"(A,A)")  "#     Zone name: ",                 zone_name
      write(*,"(A,I9)") "#     Nodes: ",                     mesh_info(1)
      write(*,"(A,I9)") "#     Cells: ",                     mesh_info(2)
      write(*,"(A,I9)") "#     Boundary nodes(if sorted): ", mesh_info(3)

      n_nodes = n_nodes + mesh_info(1)
      n_cells = n_cells + mesh_info(2)
      n_bound = n_bound + mesh_info(3)

      ! Get type of zone (structured or unstructured) 
      call Cg_Zone_Type_F(file_idx,  & ! CGNS file index number
                          base_idx,  & ! Base index number
                          zone_idx,  & ! Zone index number
                          zone_type, & ! Structured or Unstructured
                          ier)         ! Error status
        if (ier .ne. 0) then
          write(*,"(A)") "#     Failed to get zone type"
          call Cg_Error_Exit_f()
        endif

      write(*,"(A,A)") "#     Zone type is ", ZoneTypeName(zone_type)

      if (zone_type .eq. Structured) then
        write(*,"(A)") "#     Structured cgns meshed are unsupported"
        stop
      endif
      index_dim = 1

      !---------------------------!
      !   Read coordinates block  !
      !---------------------------!

      write(*,"(A,I9)") "#     Index dimension = ", index_dim

      ! Get number of coordinate arrays
      call Cg_Ncoords_F(file_idx, & ! CGNS file index number
                        base_idx, & ! Base index number
                        zone_idx, & ! Zone index number
                        n_coords, & ! Number of coordinate arrays for zone
                        ier)        ! Error status
        if (ier .ne. 0) then
          write(*,"(A)") "#       FAILED to get number of coordinate arrays"
          call Cg_Error_Exit_f()
        endif
      write(*,"(A,I9)")"#       Number of coordinate arrays for zone:", n_coords

      ! Read x, y and z coordinates
      coordinates: do coord_idx = 1, n_coords
        write(*,"(A)") "# ---------------------------"

        ! Get info about a coordinate array 
        call Cg_Coord_Info_F(file_idx,   & ! CGNS file index number
                             base_idx,   & ! Base index number
                             zone_idx,   & ! Zone index number
                             coord_idx,  & ! Coordinate array index number
                             data_type,  & ! RealSingle or RealDouble
                             coord_name, & ! Name of the coordinate array
                             ier)          ! Error status

        write(*,"(A,I3)") "#       Coord. idx:", coord_idx
        write(*,"(A,A)")  "#       Data_type: ", DataTypeName(data_type)
        write(*,"(A,A)")  "#       Name: ",      coord_name

        ! data length to be read [multizone compatible?]
        first_node = 1
        last_node  = n_nodes

        if (data_type .eq. RealSingle) then

          if(.not. allocated(data_single)) &
            allocate(data_single(first_node:last_node))
          ! Read grid coordinates (RealSingle)
          call Cg_Coord_Read_F(file_idx,    & ! CGNS file index number
                               base_idx,    & ! Base index number
                               zone_idx,    & ! Zone index number
                               coord_name,  & ! Name of the coordinate array
                               RealSingle,  & ! RealSingle or RealDouble
                               first_node,  & ! Lower range index
                               last_node,   & ! Upper range index
                               data_single, & ! Array of coordinate values
                               ier)           ! Error status
          if (ier.ne.0) then
            write(*,"(A)")"#       FAILED to read SingleReal Coord"
            call Cg_Error_Exit_f()
          endif

        elseif (data_type .eq. RealDouble) then
          if(.not. allocated(data_double)) &
            allocate(data_double(first_node:last_node))
          ! Read grid coordinates (RealDouble)
          call Cg_Coord_Read_F(file_idx,    & ! CGNS file index number
                               base_idx,    & ! Base index number
                               zone_idx,    & ! Zone index number
                               coord_name,  & ! Name of the coordinate array
                               RealDouble,  & ! RealSingle or RealDouble
                               first_node,  & ! Lower range index
                               last_node,   & ! Upper range index
                               data_double, & ! Array of coordinate values
                               ier)           ! Error status
          if (ier.ne.0) then
            write(*,"(A)")"#       FAILED to read DoubleReal Coord"
            call Cg_Error_Exit_f()
          endif

        end if

        ! append data_double to buffer array

      end do coordinates

      !---------------------!
      !   Read cells block  !
      !---------------------!

      ! Get number of element sections 
      call Cg_Nsections_F(file_idx, & ! CGNS file index number
                          base_idx, & ! Base index number
                          zone_idx, & ! Zone index number
                          n_sect,   & ! Number of element sections
                          ier)        ! Error status
      write(*,"(A)") "# ---------------------------"
      write(*,"(A,I9)") "#       Number of sections: ", n_sect

      ! Browse through all sections to read elements
      elements_connection: do sect_idx = 1, n_sect
        

        ! Get info for an element section
        ! recieves name_sect, first_cell: last_cell, n_bnd and cell_type
        call Cg_Section_Read_F(file_idx,     & ! CGNS file index number
                               base_idx,     & ! Base index number
                               zone_idx,     & ! Zone index number
                               sect_idx,     & ! Element section index
                               name_sect,    & ! Name of the Elements_t node
                               itype,        & ! Type of element
                               first_cell,   & ! Index of first element
                               last_cell,    & ! Index of last element
                               n_bnd,        & ! Index of last boundary element
                               iparent_flag, & ! if the parent data are defined
                               ier)            ! Error status

        write(*,"(A,A)")  "#         section name: ", name_sect
        write(*,"(A,A)")  "#         section type: ", ElementTypeName(itype)
        write(*,"(A,I9)") "#         first idx:",     first_cell
        write(*,"(A,I9)") "#         last idx:",      last_cell             

        ! Read element data 
        call Cg_Elements_Read_F(file_idx,       & ! CGNS file index number
                                base_idx,       & ! Base index number
                                zone_idx,       & ! Zone index number
                                sect_idx,       & ! Element section index
                                grid % cells_n, & ! Element connectivity data
                                iparent_data,   & ! For boundary or interface 
                                ier)              ! Error status

        write(*,"(A)") "# ---------------------------"

      end do elements_connection

    end do zones

  end do bases

  grid % n_nodes     = n_nodes
  grid % n_cells     = n_cells
  grid % n_bnd_cells = n_bound ! ???

  ! Allocate memory =--> carefull, there is no checking!
    ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)


  !do c = 1, grid % n_cells
  !  write(*,"(A)") "# Cell: ", c
  !  write(*,"(A)") "# Nodes: ", (grid % cells_n(n,c), n=1,8)
  !end do

stop

  end subroutine
