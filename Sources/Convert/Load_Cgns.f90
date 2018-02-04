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
  use Grid_Mod
  use Tokenizer_Mod
  use Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in
  integer*8         :: c, i, bc, base, block, sect, coord
!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = ".cgns"

  file_name = name_in

  !---------------------------------------!
  !                                       !
  !   First run: just read info from DB   !
  !                                       !
  !---------------------------------------!
  call Initialize_Counters

  ! Open a CGNS file (->file_id)
  call Cgns_Mod_Open_File

  ! Read number of CGNS bases in file_id (->n_bases)
  call Cgns_Mod_Read_Number_Of_Bases_In_File

  !------------------------------!
  !   Browse through all bases   !
  !------------------------------!
  do base = 1, n_bases

    ! Read CGNS base information in base_id
    call Cgns_Mod_Read_Base_Info(base)

    !-------------------------------------------!
    !   Browse through all blocks in the base   !
    !-------------------------------------------!
    call Cgns_Mod_Read_Number_Of_Blocks_In_Base(base)
    do block = 1, cgns_base(base) % n_blocks

      ! Read block_id information (-> n_nodes, n_cells)
      call Cgns_Mod_Read_Block_Info(base, block)

      ! Read type of block_id (->structured or unstructured)
      call Cgns_Mod_Read_Block_Type(base, block)

      !--------------------------------------------!
      !   Browse through all boundary conditions   !
      !--------------------------------------------!
      call Cgns_Mod_Read_Number_Of_Bnd_Conds_In_Block(base, block)
      do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
        call Cgns_Mod_Read_Bnd_Conds_Info(base, block, bc)
      end do

      !-----------------------------------------!
      !   Browse through all element sections   !
      !-----------------------------------------!
      call Cgns_Mod_Read_Number_Of_Element_Sections(base, block)
      do sect = 1, cgns_base(base) % block(block) % n_sects

        ! Read info for an element section (-> cell_type, first_cell, last_cell)
        call Cgns_Mod_Read_Section_Info(base, block, sect)

      end do ! elements sections

    end do ! blocks
  end do ! bases

  !------------------------!
  !   First run: results   !
  !------------------------!

  print *, "# First run finished!"
  print *, "# Found nodes: ",           cnt_nodes
  print *, "# Found cells: ",           cnt_cells
  print *, "# Found hex cells: ",       cnt_hex
  print *, "# Found pyramids cells: ",  cnt_pyr
  print *, "# Found prism cells: ",     cnt_wed
  print *, "# Found tetra cells: ",     cnt_tet
  print *, "# Found triangles faces: ", cnt_tri 
  print *, "# Found quads faces: ",     cnt_qua
  if (cnt_qua + cnt_tri .eq. 0) then
    print *, "# No boundary faces were found !"
    stop
  end if
  print *, "# Found bounary conditions faces: ", cnt_qua + cnt_tri

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  grid % n_nodes     = cnt_nodes
  grid % n_cells     = cnt_cells
  grid % n_bnd_cells = cnt_tri + cnt_qua

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)

  !-------------------------------!
  !   Allocate memory for arrays  !
  !-------------------------------!

  allocate( cgns_hex_cell_n(1:8, 1:cnt_hex) )
  allocate( cgns_pyr_cell_n(1:5, 1:cnt_pyr) )
  allocate( cgns_wed_cell_n(1:6, 1:cnt_wed) )
  allocate( cgns_tet_cell_n(1:4, 1:cnt_tet) )
  allocate( cgns_qua_face_n(1:4, 1:cnt_qua) )
  allocate( cgns_tri_face_n(1:3, 1:cnt_tri) )

  !-------------------------------------!
  !                                     !
  !   Second run: read arrays from DB   !
  !                                     !
  !-------------------------------------!
  call Initialize_Counters

  print *, "# Filling arrays.."

  !------------------------------!
  !   Browse through all bases   !
  !------------------------------!
  do base = 1, n_bases

    !-------------------------------!
    !   Browse through all blocks   !
    !-------------------------------!
    do block = 1, cgns_base(base) % n_blocks

      !---------------------------!
      !   Read coordinates block  !
      !---------------------------!

      ! Reads number of coordinates arrays from block_id (->n_coords)
      call Cgns_Mod_Read_Number_Of_Coordinates_In_Block(base, block)

      ! Read x, y and z coordinates
      do coord = 1, cgns_base(base) % block (block) % n_coords

        ! Reads coord_name of coord_id(->coord_name)
        call Cgns_Mod_Read_Coordinate_Info(base, block, coord)

        ! Read grid coordinates (-> x_, y_, z_coord)
        call Cgns_Mod_Read_Coordinate_Array(base, block, coord, grid)

      end do ! coordinates

      !---------------------!
      !   Read cells block  !
      !---------------------!

      ! Browse through all sections to read elements
      do sect = 1, cgns_base(base) % block(block) % n_sects

        ! Read element data (count HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
        call Cgns_Mod_Read_Section_Connections(base, block, sect)

      end do ! elements sections

    end do ! blocks

  end do ! bases

stop

  !--------------------------------------!
  !  Conversion to T-FlowS B.C. format   !
  !--------------------------------------!

  ! Reorder elements connectivity according to neu_mod

  ! HEXA_8 (-> f8n)
  do c = lbound(cgns_hex_cell_n,dim=2), &
    ubound(cgns_hex_cell_n,dim=2) - lbound(cgns_hex_cell_n,dim=2) + 1

    i = cgns_hex_cell_n(4, c)
    cgns_hex_cell_n(4, c) = cgns_hex_cell_n(3, c)
    cgns_hex_cell_n(3, c) = i
    i = cgns_hex_cell_n(8, c)
    cgns_hex_cell_n(8, c) = cgns_hex_cell_n(7, c)
    cgns_hex_cell_n(7, c) = i
  end do

  ! PYRA_5 (-> f5n)
  do c =  lbound(cgns_pyr_cell_n,dim=2), &
    ubound(cgns_pyr_cell_n,dim=2) - lbound(cgns_pyr_cell_n,dim=2) + 1

    i = cgns_pyr_cell_n(4, c)
    cgns_pyr_cell_n(4, c) = cgns_pyr_cell_n(3, c)
    cgns_pyr_cell_n(3, c) = i
  end do

  ! PENTA_6 (-> f6n)
  ! no changes

  ! TETRA_4 (-> f4n)
  ! no changes

stop

end subroutine
