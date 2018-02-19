!==============================================================================!
  subroutine Cgns_Mod_Write_Field_Par(base, block, solution, field, grid, &
    input_array, input_name)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field_id  [parallel vesion]     !
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!   Connections:  |-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|...!
!----------------------------------[Modules]-----------------------------------!
  use par_mod, only: this_proc
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base, block, solution, field
  type(Grid_Type)      :: grid
  real                 :: input_array(1:grid % n_cells)
  character(len=*)     :: input_name
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id        ! base index number
  integer              :: block_id       ! block index number
  integer              :: solution_id    ! solution index
  integer              :: field_id       ! field index
  character(len=80)    :: field_name     ! name of the FlowSolution_t node
  integer              :: cell_type
  integer              :: cnt            ! cells of cell_type
  real, allocatable    :: field_array(:) ! field array
  integer              :: i, j, k, c
  integer              :: error
!==============================================================================!

  ! Set input parameters
  base_id       = base
  block_id      = block
  solution_id   = solution
  field_id      = field

  field_name = trim(input_name)

  !----------------------------------------------------!
  !   Add an empty field node to FlowSolution_t node   !
  !----------------------------------------------------!

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

  !------------------------------------------------!
  !   Mapping 1:NC -> Connection structure above   !
  !------------------------------------------------!

  ! Find first and last cells of cell_type
  do cell_type = 1, 4
    ! cells of cell_type
    if (cell_type.eq.1) cnt = cnt_hex
    if (cell_type.eq.2) cnt = cnt_wed
    if (cell_type.eq.3) cnt = cnt_pyr
    if (cell_type.eq.4) cnt = cnt_tet

    i = cnt ! cnt_hex/pyr/wed/tet on this_proc
    call Cgns_Mod_Get_Arrays_Dimensions_Par(j, i)

    i = j + cnt_cells
    j = i + cnt - 1

    if (cnt.ne.0) then

      allocate(field_array(i:j), stat = error); field_array = 0

      if (error .ne. 0) then
         print*, '*FAILED* to allocate ', "field_array"
         call Cg_Error_Exit_F()
      endif

      ! copy input array to field_array
      k = i
      do c = 1, grid % n_cells
        if     (cell_type.eq.1 .and. grid % cells_n_nodes(c).eq.8) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (cell_type.eq.2 .and. grid % cells_n_nodes(c).eq.6) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (cell_type.eq.3 .and. grid % cells_n_nodes(c).eq.5) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (cell_type.eq.4 .and. grid % cells_n_nodes(c).eq.4) then
          field_array(k) = input_array(c)
          k = k + 1
        end if
      end do

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

      call wait
      deallocate(field_array)

      c = cnt
      call wait
      call IglSum(c)
      cnt_cells = cnt_cells + c

      ! Print some info
      if(verbose .and. this_proc.eq.1) then
        print *, '#           ---------------------------------'
        print *, '#           Field name: ',  field_name
        print *, '#           Field idx:    ', field_id
        print *, '#           ---------------------------------'
      end if
      if(verbose) then
        print *, '#           First cell:', i, " (P:",this_proc,")"
        print *, '#           Last cell: ', j, " (P:",this_proc,")"
      end if

    end if

  end do

  cnt_cells = 0

  end subroutine