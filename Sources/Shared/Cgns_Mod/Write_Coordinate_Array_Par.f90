!==============================================================================!
  subroutine Cgns_Mod_Write_Coordinate_Array_Par(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Gets n_sects from block
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  ! 
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, coord
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id         ! base index number
  integer           :: block_id        ! block index number
  integer           :: coord_id        ! coord. array index number
  character(len=80) :: coord_name      
  integer           :: i               ! lower range index
  integer           :: j               ! upper range index
  real, allocatable :: coordinates(:)  ! array of coordinate values
  integer           :: error           ! error status
  integer           :: first_node      ! look at array structure at the header
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  coord_id   = coord
  coord_name = cgns_base(base) % block(block) % coord_name(coord)

  !---------- fetch x, y, z arrays structure
  i = grid % n_nodes
  call Cgns_Mod_Get_Arrays_Dimensions_Par(first_node, i)

  i = first_node
  j = first_node + grid % n_nodes - 1

  allocate(coordinates(i:j), stat = error); coordinates = 0
  if (error .ne. 0) then
     print*, '*FAILED* to allocate ', "coordinates"
     call cgp_error_exit_f()
  endif


  ! Fetch received parameters
  select case (coord_id)
    case (1)
      coordinates(i:j) = grid % xn(:)
    case (2)
      coordinates(i:j) = grid % yn(:)
    case (3)
      coordinates(i:j) = grid % zn(:)
  end select

  !---------- Create empty coord_name node in DB
  call Cgp_Coord_Write_F( &
    file_id,              &
    base_id,              &
    block_id,             &
    RealDouble,           &
    coord_name,           &
    coord_id,             &
    error)

    if (error.ne.CG_OK) then
      print*,'*FAILED* to create empty: ', trim(coord_name)
      call cgp_error_exit_f()
    endif
  !---------- Fill that node with grid coordinates
  call Cgp_Coord_Write_Data_F( &
     file_id,                  &
     base_id,                  &
     block_id,                 &
     coord_id,                 &
     i,                        &
     j,                        &
     coordinates(i:j),         &
     error)

    if (error .ne. 0) then
           print *, "# Failed to fill: ", trim(coord_name)
       call cgp_error_exit_f()
    endif
  deallocate(coordinates)

  ! Print some info
  if(verbose .and. this_proc.eq.1) then
    print *, '#         Coord array: ', coord_name
  end if
  if(verbose.and.coord_id.eq.1) then
    print *, '#         Number of nodes: ', j - i + 1, " (P:",this_proc,")"
    print *, '#         First node:', i, " (P:",this_proc,")"
    print *, '#         Last node: ', j, " (P:",this_proc,")"
  end if

  end subroutine