!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Zones_In_Base
!------------------------------------------------------------------------------!
!   Gets n_zones from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

 ! Get number of zones in base
  call Cg_Nzones_F(file_id, & ! cgns file index number
                   base_id, & ! base index number
                   n_zones, & ! number of zones present in base
                   ier)       ! error status
  if (ier .ne. 0) then
    print *, "#   FAILED to get zones number"
    call Cg_Error_Exit_F()
  endif

  print *, "#   Total zones:", n_zones

  end subroutine
