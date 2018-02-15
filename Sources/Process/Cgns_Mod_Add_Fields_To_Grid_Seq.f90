!==============================================================================!
  subroutine Cgns_Mod_Add_Fields_To_Grid_Seq(grid, name_save)
!------------------------------------------------------------------------------!
!   Closes file_id file [sequential vesion]                                    !
!------------------------------------------------------------------------------!
  use all_mod
  use Grid_Mod
  use Cgns_Mod
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: store_name
  integer           :: c, base, block, sect, coord, mode, solution, field
!==============================================================================!

  print *, "# subroutine Cgns_Mod_Save_Grid_Seq"

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  !--------------------------!
  !   Open file for modify   !
  !--------------------------!
  call Name_File(0, file_name, '.cgns')

  mode = CG_MODE_MODIFY
  call Cgns_Mod_Open_File_Seq(mode)

  !-----------------!
  !                 !
  !   Bases block   !
  !                 !
  !-----------------!
  n_bases = 1
  allocate(cgns_base(n_bases))

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

  !--------------------!
  !                    !
  !   Solution block   !
  !                    !
  !--------------------!

  cgns_base(base) % block(block) % n_solutions = 1

  allocate(cgns_base(base) % block(block) % solution( &
    cgns_base(base) % block(block) % n_solutions))
  solution = 1

  cgns_base(base) % block(block) % solution(solution) % name = 'FlowSolution'

  call Cgns_Mod_Write_Solution_Info_Seq(base, block, solution, grid)

  !-----------------!
  !                 !
  !   Field block   !
  !                 !
  !-----------------!
  cgns_base(base) % block(block) % solution(solution) % n_fields = 1
  allocate(cgns_base(base) % block(block) % solution(solution) % field( &
    cgns_base(base) % block(block) % solution(solution) % n_fields))

  field = 1


  ! Close file
  call Cgns_Mod_Close_File_Seq

  end subroutine
