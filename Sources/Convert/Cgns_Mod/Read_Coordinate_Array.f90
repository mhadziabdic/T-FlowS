!==============================================================================!
  subroutine Cgns_Mod_Read_Coordinate_Array(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Read grid coordinates (RealDouble)                                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8       :: base, block, coord
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer*8         :: base_id
  integer*8         :: block_id
  character(len=80) :: coord_name
  integer*8         :: error
  integer*8         :: i, j
!==============================================================================!
!   Description of arguments for the called CGNS function:
!
!   Cg_Coord_Read_F(file_id,        &  ! cgns file index number
!                   base,           &  ! base index number
!                   block,          &  ! block index number
!                   coord_name,     &  ! name of the coordinate array
!                   RealDouble,     &  ! realsingle or realdouble
!                   i,              &  ! lower range index
!                   j,              &  ! upper range index
!                   buffer_double,  &  ! array of coordinate values
!                   error)               ! error status
!------------------------------------------------------------------------------!

  ! Set input parameters
  base_id  = base
  block_id = block

  i = 1
  j = cgns_base(base) % block(block) % mesh_info(1)

  allocate(buffer_double(i:j))

  ! Read grid x coordinates
  call Cg_Coord_Read_F(file_id,        &
                       base,           &
                       block,          &
                       coord_name,     &
                       RealDouble,     &
                       i,              &
                       j,              &
                       buffer_double,  &
                       error)

  if (error.ne.0) then
    print *, "# Failed to read DoubleReal Coord", coord_name
    call Cg_Error_Exit_F()
  endif

  select case (coord)
    case (1)
      i = cnt_x + 1
      j = cnt_x + j
      grid % xn(i:j) = buffer_double(:)
      cnt_x = j
    case (2)
      i = cnt_y + 1
      j = cnt_y + j
      grid % yn(i:j) = buffer_double(:)
      cnt_y = j
    case (3)
      i = cnt_z + 1
      j = cnt_z + j
      grid % zn(i:j) = buffer_double(:)
      cnt_z = j
  end select

  deallocate(buffer_double)

  end subroutine
