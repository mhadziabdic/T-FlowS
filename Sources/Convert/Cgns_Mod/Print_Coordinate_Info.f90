!==============================================================================!
  subroutine Cgns_Mod_Print_Coordinate_Info(base, block)
!------------------------------------------------------------------------------!
!   Reads coord_name of coord_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer*8 :: ier
!==============================================================================!
!   Description of arguments for the CGNS function call:
!
!   Cg_Coord_Info_F(file_id,          & 
!                   base,             & 
!                   block,            & 
!                   coord_id,         & 
!                   coord_data_type,  & 
!                   coord_name,       & 
!                   ier)                
!------------------------------------------------------------------------------!

  !Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,          &
                       base,             &
                       block,            &
                       coord_id,         &
                       coord_data_type,  &
                       coord_name,       &
                       ier)             
  if (ier .ne. 0) then
    print *, "# Failed to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  print *, "# Coord. id:",         coord_id
  print *, "# Coord. Data Type: ", DataTypeName(coord_data_type)
  print *, "# Name: ",             coord_name

  end subroutine
