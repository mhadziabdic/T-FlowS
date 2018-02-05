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
  integer*8         :: c, cs, ce, i, j, bc, base, block, sect, coord
  integer*8         :: cgns_1, cgns_2, cgns_3, cgns_4, cgns_5, cell_type
!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = '.cgns'

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

  print *, '# First run finished!'
  print *, '# - number of nodes: ',           cnt_nodes
  print *, '# - number of cells: ',           cnt_cells
  print *, '# - number of hex cells: ',       cnt_hex
  print *, '# - number of pyramids cells: ',  cnt_pyr
  print *, '# - number of prism cells: ',     cnt_wed
  print *, '# - number of tetra cells: ',     cnt_tet
  print *, '# - number of triangles faces: ', cnt_tri 
  print *, '# - number of quads faces: ',     cnt_qua
  if (cnt_qua + cnt_tri .eq. 0) then
    print *, '# No boundary faces were found !'
    stop
  end if
  print *, '# - number of bounary conditions faces: ', cnt_qua + cnt_tri
  print *, '# - number of bounary conditions: ', cnt_bnd_conds
  do i = 1, cnt_bnd_conds
    print *, bnd_cond_names(i)
  end do 

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  grid % n_nodes     = cnt_nodes
  grid % n_cells     = cnt_cells
  grid % n_bnd_cells = cnt_tri + cnt_qua

  call Allocate_Memory(grid)

  !-------------------------------------!
  !                                     !
  !   Second run: read arrays from DB   !
  !                                     !
  !-------------------------------------!
  call Initialize_Counters

  print *, '# Filling arrays..'

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
      cnt_block_bnd_cells = 0

      ! Browse through all sections to read elements
      do sect = 1, cgns_base(base) % block(block) % n_sects

        ! Read element data (count HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
        call Cgns_Mod_Read_Section_Connections(base, block, sect, grid)

      end do ! elements sections

      cnt_nodes = cnt_nodes + cgns_base(base) % block(block) % mesh_info(1)
      cnt_cells = cnt_cells + cgns_base(base) % block(block) % mesh_info(2)
      cnt_bnd_cells = cnt_bnd_cells + cnt_block_bnd_cells

    end do ! blocks
  end do ! bases

  print *, '# Grid metrix after reading:'
  print *, '# - number of nodes:         ', cnt_nodes
  print *, '# - number of cells:         ', cnt_cells   
  print *, '# - number of boundary cells:', cnt_bnd_cells

  call Merge_Nodes(grid)

  !-----------------------------------------------------------------!
  !   Correct boundary conditions directions for hexahedral cells   !
  !   (They are not the same in CGNS and Gambit's neutral format.)  !
  !-----------------------------------------------------------------!
  cnt_bnd_cells = 0
  do base = 1, n_bases
    do block = 1, cgns_base(base) % n_blocks
      do sect = 1, cgns_base(base) % block(block) % n_sects
        cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type

        if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then
          cs = cgns_base(base) % block(block) % section(sect) % first_cell
          ce = cgns_base(base) % block(block) % section(sect) % last_cell

          do c = cs, ce
            cgns_1 = grid % cells_bnd_color(1,c)
            cgns_2 = grid % cells_bnd_color(2,c)
            cgns_3 = grid % cells_bnd_color(3,c)
            cgns_4 = grid % cells_bnd_color(4,c)
            cgns_5 = grid % cells_bnd_color(5,c)
            grid % cells_bnd_color(4,c) = cgns_5
            grid % cells_bnd_color(3,c) = cgns_4
            grid % cells_bnd_color(2,c) = cgns_3
            grid % cells_bnd_color(1,c) = cgns_2
            grid % cells_bnd_color(5,c) = cgns_1

            do j = 1, 6
              if( grid % cells_bnd_color(j,c) .ne. 0 ) then 
                cnt_bnd_cells = cnt_bnd_cells + 1
              end if
            end do 
          end do 

        end if

      end do ! elements sections
    end do ! blocks
  end do ! bases
  print *, '# - number of boundary cells:', cnt_bnd_cells

  end subroutine
