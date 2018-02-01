!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Bases_In_File
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get number of CGNS bases in file
  call Cg_Nbases_F(file_id, & ! cgns file index number
                   n_bases, & ! number of bases present in the CGNS file
                   ier)       ! error status

    if (ier .ne. 0) then
      print *, "# Failed to get bases number"
      call Cg_Error_Exit_F()
    endif

  print *, "# Number of bases: ", n_bases

  end subroutine
