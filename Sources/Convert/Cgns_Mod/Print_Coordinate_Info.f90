!==============================================================================!
  subroutine Cgns_Mod_Print_Coordinate_Info
!------------------------------------------------------------------------------!
!   Reads coord_name of coord_id
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,         & ! cgns file index number
                       base_id,         & ! base index number
                       zone_id,         & ! zone index number
                       coord_id,        & ! Coordinate array index number
                       coord_data_type, & ! realsingle or realdouble
                       coord_name,      & ! name of the coordinate array
                       ier)               ! error status

  if (ier .ne. 0) then
    print *, "# Failed to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  print *, "# Coord. id:",         coord_id
  print *, "# Coord. Data Type: ", DataTypeName(coord_data_type)
  print *, "# Name: ",             coord_name

  end subroutine
