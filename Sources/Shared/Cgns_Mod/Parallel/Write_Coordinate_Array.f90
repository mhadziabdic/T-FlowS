!==============================================================================!
  subroutine Write_Coordinate_Array(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Writes grid coordinates (RealDouble) [parallel vesion]                     !
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
  coord_name = cgns_base(base_id) % block(block_id) % coord_name(coord_id)

  !----------------------------------------------!
  !   Create empty coord_name node in DB block   !
  !----------------------------------------------!

  ! Fetch coordinates array dimensions
  i = grid % n_nodes
  call Get_Arrays_Dimensions(first_node, i)

  i = first_node
  j = first_node + grid % n_nodes - 1

  allocate(coordinates(i:j), stat = error); coordinates = 0
  if (error .ne. 0) then
     print*, '*FAILED* to allocate ', "coordinates"
     call Cgp_Error_Exit_F()
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

  call Cgp_Coord_Write_F( & !(in )
    file_id,              & !(in )
    base_id,              & !(in )
    block_id,             & !(in )
    RealDouble,           & !(in )
    coord_name,           & !(in )
    coord_id,             & !(out)
    error)                  !(out)

  if (error.ne.CG_OK) then
    print*,'*FAILED* to create empty: ', trim(coord_name)
    call Cgp_Error_Exit_F()
  endif

  !------------------------------------------!
  !   Fill empty coord_id node in DB block   !
  !------------------------------------------!

  ! Fill coord_id node with grid coordinates
  call Cgp_Coord_Write_Data_F( & !(in )
     file_id,                  & !(in )
     base_id,                  & !(in )
     block_id,                 & !(in )
     coord_id,                 & !(in )
     i,                        & !(in )
     j,                        & !(in )
     coordinates(i:j),         & !(in )
     error)                      !(out)

  if (error .ne. 0) then
         print *, "# Failed to fill: ", trim(coord_name)
     call Cgp_Error_Exit_F()
  endif

  deallocate(coordinates)

  ! Print some info
  if(verbose .and. this_proc.eq.1) then
    print *, '#         Coord array: ', coord_name
  end if
  if(verbose.and.coord_id.eq.1) then
    print *, '#         Number of nodes: ', j - i + 1, " (P:",this_proc,")"
    print *, '#         First node:', i,               " (P:",this_proc,")"
    print *, '#         Last node: ', j,               " (P:",this_proc,")"
  end if

  end subroutine