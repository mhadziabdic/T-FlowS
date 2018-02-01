!==============================================================================!
  subroutine Cgns_Mod_Read_Zone_Type
!------------------------------------------------------------------------------!
!   Prints zone type                                                           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Read type of zone
  call Cg_Zone_Type_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      zone_id,   & ! zone index number
                      zone_type, & ! structured or unstructured
                      ier)         ! error status
  if (ier .ne. 0) then
    print *, "# Failed to get zone type"
    call Cg_Error_Exit_F()
  endif

  print *, "# Zone type is ", ZoneTypeName(zone_type)

  if (zone_type .eq. Structured) then
    print *, "# Structured cgns meshed are unsupported"
    stop
  endif

  print *, "# Index dimension = ", 1

  end subroutine
