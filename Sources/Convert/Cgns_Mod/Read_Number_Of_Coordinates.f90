!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Coordinates
!------------------------------------------------------------------------------!
!   Reads number of coordinates arrays from zone_id
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get number of coordinate arrays (1 for unstructure, 3 for structured)
  call Cg_Ncoords_F(file_id,  & ! cgns file index number
                    base_id,  & ! base index number
                    zone_id,  & ! zone index number
                    n_coords, & ! number of coordinate arrays for zone
                    ier)        ! error status

  if (ier .ne. 0) then
    print *, "#       FAILED to get number of coordinate arrays"
    call Cg_Error_Exit_F()
  endif

  print *, "#       Number of coordinate arrays for zone:", n_coords

  end subroutine
