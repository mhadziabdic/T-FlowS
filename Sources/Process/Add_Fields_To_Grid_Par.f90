!==============================================================================!
  subroutine Add_Fields_To_Grid_Par(grid, name_save)
!------------------------------------------------------------------------------!
!   Closes file_id file [Parallel vesion]                                      !
!------------------------------------------------------------------------------!
  use all_mod
  use par_mod, only: this_proc
  use pro_mod
  use Grid_Mod
  use Cgns_Mod
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: store_name
  integer           :: base
  integer           :: block
  integer           :: solution
  integer           :: field
!==============================================================================!

  print *, "# subroutine Cgns_Mod_Save_Grid_Par"

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  !--------------------------!
  !   Open file for modify   !
  !--------------------------!
  call Name_File(0, file_name, '.cgns')

  file_mode = CG_MODE_MODIFY
  call Cgns_Mod_Open_File_Par(file_mode)

  call Initialize_Counters

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
  cgns_base(base) % block(block) % solution(solution) % sol_type = CellCenter

  call Cgns_Mod_Write_Solution_Info_Par(base, block, solution)

  !-----------------!
  !                 !
  !   Field block   !
  !                 !
  !-----------------!
  cgns_base(base) % block(block) % solution(solution) % n_fields = 4
  allocate(cgns_base(base) % block(block) % solution(solution) % field( &
    cgns_base(base) % block(block) % solution(solution) % n_fields))

  !-------------------!
  !   Scalar fields   !
  !-------------------!

  field = 1
  cgns_base(base)%block(block)%solution(solution)% &
    field(field)%name = 'Pressure'
  call Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    P % n(1:grid % n_cells))

  !-------------------!
  !   Vector fields   !
  !-------------------!

  field = 2
  cgns_base(base)%block(block)%solution(solution)% &
    field(field)%name = 'VelocityX'
  call Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    U % n(1:grid % n_cells))

  field = 3
  cgns_base(base)%block(block)%solution(solution)% &
    field(field)%name = 'VelocityY'
  call Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    V % n(1:grid % n_cells))

  field = 4
  cgns_base(base)%block(block)%solution(solution)% &
    field(field)%name = 'VelocityZ'
  call Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    W % n(1:grid % n_cells))

  !----------------------!
  !   Pack in function   !
  !----------------------!

!!   go to base node
!      call cg_goto_f(index_file,index_base,ier,'end')
!!   write descriptor node (user can give any name)
!      text1='Supersonic vehicle with landing gear'
!      text2='M=4.6, Re=6 million'
!      textstring=text1//char(10)//text2
!      call cg_descriptor_write_f('Information',textstring,ier)

  ! Close DB
  call Cgns_Mod_Close_File_Par

  if (this_proc.eq.1) &
    print *, 'Successfully added fields to ', trim(problem_name)

  deallocate(cgns_base)

  ! Restore the name
  problem_name = store_name

  end subroutine
