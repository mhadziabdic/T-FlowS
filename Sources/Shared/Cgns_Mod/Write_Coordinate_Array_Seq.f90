!==============================================================================!
  subroutine Cgns_Mod_Write_Coordinate_Array_Seq(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Gets n_sects from block
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, coord
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id         ! base index number
  integer           :: block_id        ! block index number
  character(len=80) :: coord_name      
  integer           :: i               ! lower range index
  integer           :: j               ! upper range index
  real, allocatable :: coordinates(:)  ! array of coordinate values
  integer           :: error           ! error status
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  coord_name = cgns_base(base) % block(block) % coord_name(coord)

  i = 1
  j = cgns_base(base) % block(block) % mesh_info(1)
  allocate(coordinates(i:j))

  ! Fetch received parameters
  select case (coord)
    case (1)
      coordinates(i:j) = grid % xn(i:j)
    case (2)
      coordinates(i:j) = grid % yn(i:j)
    case (3)
      coordinates(i:j) = grid % zn(i:j)
  end select

  !---------- Write grid coordinates 
  call Cg_Coord_Write_F( &
    file_id,             &
    base_id,             &
    block_id,            &
    RealDouble,          &
    coord_name,          &
    coordinates(i:j),    &
    coord,               &
    error)

    if (error .ne. 0) then
           print *, "# Failed to write: ", trim(coord_name)
       call cg_error_exit_f()
    endif
  deallocate(coordinates)

  end subroutine