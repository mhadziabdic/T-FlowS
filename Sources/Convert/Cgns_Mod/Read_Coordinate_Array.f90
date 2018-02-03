!==============================================================================!
  subroutine Cgns_Mod_Read_Coordinate_Array(base, block)
!------------------------------------------------------------------------------!
!   Read grid coordinates (RealDouble)                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer*8 :: i, j, ier
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
!                   ier)               ! error status
!------------------------------------------------------------------------------!

  i = 1
  j = cgns_block(block) % mesh_info(1)

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
                       ier)

  if (ier.ne.0) then
    print *, "# Failed to read DoubleReal Coord", coord_name
    call Cg_Error_Exit_F()
  endif

  select case (coord_id)
    case (1)
      i = last_x + 1
      j = last_x + j
      x_coord(i:j) = buffer_double(:)
      last_x = j
    case (2)
      i = last_y + 1
      j = last_y + j
      y_coord(i:j) = buffer_double(:)
      last_y = j
    case (3)
      i = last_z + 1
      j = last_z + j
      z_coord(i:j) = buffer_double(:)
      last_z = j
  end select

  deallocate(buffer_double)

  end subroutine
