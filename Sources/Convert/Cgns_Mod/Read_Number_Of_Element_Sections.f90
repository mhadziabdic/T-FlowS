!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Element_Sections
!------------------------------------------------------------------------------!
!   Gets n_sects from zone
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get number of element sections
  call Cg_Nsections_F(file_id,  & ! cgns file index number
                      base_id,  & ! base index number
                      zone_id,  & ! zone index number
                      n_sects,  & ! number of element sections
                      ier)        ! error status

  if (ier.ne.0) then
    print *, "# Failed to read number of elements"
    call Cg_Error_Exit_F()
  endif

  print *, "# ---------------------------"
  print *, "# Number of sections: ", n_sects

  end subroutine
