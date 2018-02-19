subroutine Cgns_Mod_Save_Grid_Seq(grid, name_save)                                        !
!   Writes in sequentialy 3-D unstructured grid to files '????????'            !
!------------------------------------------------------------------------------!
!   All conventions for this lib are here:                                     !
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html                    !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Grid_Mod
  use Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: store_name
  integer           :: c, base, block, sect, coord
!==============================================================================!

  print *, "# subroutine Cgns_Mod_Save_Grid_Seq"

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  !-------------------------!
  !   Open file for write   !
  !-------------------------!
  call Name_File(0, file_name, '.cgns')

  file_mode = CG_MODE_WRITE
  call Cgns_Mod_Open_File_Seq(file_mode)

  call Initialize_Counters

  ! Count number of 3d cell type elements
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 8) cnt_hex = cnt_hex + 1
    if(grid % cells_n_nodes(c) == 6) cnt_wed = cnt_wed + 1
    if(grid % cells_n_nodes(c) == 5) cnt_pyr = cnt_pyr + 1
    if(grid % cells_n_nodes(c) == 4) cnt_tet = cnt_tet + 1
  end do

  !-----------------!
  !                 !
  !   Bases block   !
  !                 !
  !-----------------!
  n_bases = 1
  allocate(cgns_base(n_bases))

  base = 1
  cgns_base(base) % name = "Base 1"
  cgns_base(base) % cell_dim = 3
  cgns_base(base) % phys_dim = 3

  call Cgns_Mod_Write_Base_Info_Seq(base)

  !-----------------!
  !                 !
  !   Zones block   !
  !                 !
  !-----------------!

  cgns_base(base) % n_blocks = 1
  allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

  block = 1
  cgns_base(base) % block(block) % name = "Zone 1"
  cgns_base(base) % block(block) % mesh_info(1) = grid % n_nodes
  cgns_base(base) % block(block) % mesh_info(2) = grid % n_cells
  cgns_base(base) % block(block) % mesh_info(3) = 0

  call Cgns_Mod_Write_Block_Info_Seq(base, block)

  !-----------------------!
  !                       !
  !   Coordinates block   !
  !                       !
  !-----------------------!

  coord = 1
  cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateX'
  call Cgns_Mod_Write_Coordinate_Array_Seq(base, block, coord, grid)

  coord = 2
  cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateY'
  call Cgns_Mod_Write_Coordinate_Array_Seq(base, block, coord, grid)

  coord = 3
  cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateZ'
  call Cgns_Mod_Write_Coordinate_Array_Seq(base, block, coord, grid)

  !-----------------------------!
  !                             !
  !   Cells connections block   !
  !                             !
  !-----------------------------!

  cgns_base(base) % block(block) % n_sects = 4
  allocate(cgns_base(base) % block(block) % section( &
    cgns_base(base) % block(block) % n_sects))

  sect = 1
  cgns_base(base)%block(block)%section(sect)%name = 'Hexagons'
  cgns_base(base)%block(block)%section(sect)%cell_type = HEXA_8
  cgns_base(base)%block(block)%section(sect)%first_cell = 1      ! + cnt_cells
  cgns_base(base)%block(block)%section(sect)%last_cell  = cnt_hex! + cnt_cells
  call Cgns_Mod_Write_Section_Connections_Seq(base, block, sect, grid)

  sect = 2
  cgns_base(base)%block(block)%section(sect)%name = 'Pyramids'
  cgns_base(base)%block(block)%section(sect)%cell_type = PYRA_5
  cgns_base(base)%block(block)%section(sect)%first_cell = 1     ! + cnt_cells
  cgns_base(base)%block(block)%section(sect)%last_cell = cnt_pyr! + cnt_cells
  call Cgns_Mod_Write_Section_Connections_Seq(base, block, sect, grid)

  sect = 3
  cgns_base(base)%block(block)%section(sect)%name = 'Wedges'
  cgns_base(base)%block(block)%section(sect)%cell_type = PENTA_6
  cgns_base(base)%block(block)%section(sect)%first_cell = 1     ! + cnt_cells
  cgns_base(base)%block(block)%section(sect)%last_cell = cnt_wed! + cnt_cells
  call Cgns_Mod_Write_Section_Connections_Seq(base, block, sect, grid)

  sect = 4
  cgns_base(base)%block(block)%section(sect)%name = 'Tetrahedrons'
  cgns_base(base)%block(block)%section(sect)%cell_type = TETRA_4
  cgns_base(base)%block(block)%section(sect)%first_cell = 1     ! + cnt_cells
  cgns_base(base)%block(block)%section(sect)%last_cell = cnt_tet! + cnt_cells
  call Cgns_Mod_Write_Section_Connections_Seq(base, block, sect, grid)

  ! Close DB
  call Cgns_Mod_Close_File_Seq

  print *, 'Successfully wrote unstructured grid to file ', trim(problem_name)

  deallocate(cgns_base)

  ! Restore the name
  problem_name = store_name

end subroutine