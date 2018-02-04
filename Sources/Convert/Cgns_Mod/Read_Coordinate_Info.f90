!==============================================================================!
  subroutine Cgns_Mod_Read_Coordinate_Info(base, block, coord)
!------------------------------------------------------------------------------!
!   Reads coord_name of coord_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block, coord
!-----------------------------------[Locals]-----------------------------------!
  integer*8         :: base_id          ! base index number
  integer*8         :: block_id         ! block index number
  integer*8         :: coord_id            
  integer*8         :: coord_data_type     
  character(len=80) :: coord_name
  integer*8         :: error            ! error status
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  coord_id = coord
 
  ! Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,          &
                       base_id,          &
                       block_id,         &
                       coord_id,         &
                       coord_data_type,  &
                       coord_name,       &
                       error)             
  if (error .ne. 0) then
    print *, "# Failed to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % coord_name(coord) = coord_name

  print *, "# Coord. id:",         coord_id
  print *, "# Coord. Data Type: ", DataTypeName(coord_data_type)
  print *, "# Name: ",             coord_name

  end subroutine
