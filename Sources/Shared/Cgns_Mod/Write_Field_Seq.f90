!==============================================================================!
  subroutine Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    array)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field_id  [sequential vesion] !
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base, block, solution, field
  type(Grid_Type)      :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id     ! base index number
  integer              :: block_id    ! block index number
  integer              :: solution_id ! solution index
  integer              :: field_id    ! field index
  character(len=80)    :: field_name  ! name of the FlowSolution_t node
  integer              :: error
  real                 :: array(1:grid % n_cells)
!==============================================================================!

  ! Set input parameters
  base_id       = base
  block_id      = block
  solution_id   = solution
  field_id      = field

  field_name = trim(cgns_base(base_id)%block(block_id)%solution(solution_id)% &
                    field(field_id)%name)

  ! Add field to FlowSolution_t node 
  call Cg_Field_Write_F( & !(in )
    file_id,             & !(in )
    base_id,             & !(in )
    block_id,            & !(in )
    solution_id,         & !(in )
    RealDouble,          & !(in )
    field_name,          & !(in )
    array,               & !(in )
    field_id,            & !(out)
    error)                 !(out)

  if (error .ne. 0) then
    print *, "# Failed to write field", trim(field_name)
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose ) then
    print *, '#           ---------------------------------'
    print *, '#           Field name: ',  field_name
    print *, '#           Field idx:    ', field_id
    print *, '#           ---------------------------------'
  end if


  end subroutine