!==============================================================================!
  subroutine Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    input_array)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field_id  [parallel vesion]     !
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  !
!----------------------------------[Modules]-----------------------------------!
  use par_mod, only: this_proc
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base, block, solution, field
  type(Grid_Type)      :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id        ! base index number
  integer              :: block_id       ! block index number
  integer              :: solution_id    ! solution index
  integer              :: field_id       ! field index
  character(len=80)    :: field_name     ! name of the FlowSolution_t node
  integer              :: error
  real                 :: input_array(1:grid % n_cells)
  real, allocatable    :: field_array(:) ! array of coordinate values
  integer              :: i              ! lower range index
  integer              :: j              ! upper range index
  integer              :: first_node     ! look at array structure at the header
!==============================================================================!

  ! Set input parameters
  base_id       = base
  block_id      = block
  solution_id   = solution
  field_id      = field

  field_name = trim(cgns_base(base_id)%block(block_id)%solution(solution_id)% &
                    field(field_id)%name)

  !----------------------------------------------------!
  !   Add an empty field node to FlowSolution_t node   !
  !----------------------------------------------------!

  ! Fetch coordinates field_array dimensions
  !i = grid % n_cells
  !call Cgns_Mod_Get_Arrays_Dimensions_Par(first_node, i)

  !i = first_node
  !j = first_node + grid % n_cells - 1

  i = minval(tflows_2_cgns_cells(1:grid % n_cells))
  j = maxval(tflows_2_cgns_cells(1:grid % n_cells))
  print *, "pid=", this_proc,"i=", i, "j= ", j

  allocate(field_array(i:j), stat = error); field_array = 0

  if (error .ne. 0) then
     print*, '*FAILED* to allocate ', "field_array"
     call Cgp_Error_Exit_F()
  endif

  ! copy input array to field_array
  field_array(i:j) = input_array(1:grid % n_cells)
  
  call Cgp_Field_Write_F( & !(in )
    file_id,              & !(in )
    base_id,              & !(in )
    block_id,             & !(in )
    solution_id,          & !(in )
    RealDouble,           & !(in )
    field_name,           & !(in )
    field_id,             & !(out)
    error)                  !(out)

  if (error .ne. 0) then
    print *, "# Failed to create empty ", trim(field_name)
    call Cgp_Error_Exit_F()
  endif

  !---------------------------------------!
  !   Fill empty field_name in DB block   !
  !---------------------------------------!
  call Cgp_Field_Write_Data_F( & !(in )
    file_id,                   & !(in )
    base_id,                   & !(in )
    block_id,                  & !(in )
    solution_id,               & !(in )
    field_id,                  & !(in )
    i,                         & !(in )
    j,                         & !(in )
    field_array,               & !(in )
    error)                       !(out)

  if (error .ne. 0) then
    print *, "# Failed to fill ", trim(field_name)
    call Cgp_Error_Exit_F()
  endif

  deallocate(field_array)

  ! Print some info
  if(verbose .and. this_proc.eq.1) then
    print *, '#           ---------------------------------'
    print *, '#           Field name: ',  field_name
    print *, '#           Field idx:    ', field_id
    print *, '#           ---------------------------------'
  end if


  end subroutine