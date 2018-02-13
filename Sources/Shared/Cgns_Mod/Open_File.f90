!==============================================================================!
  subroutine Cgns_Mod_Open_File(mode)
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: mode
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  print *, "# Reading the file:", trim(file_name)

  ! Open a CGNS file
  call Cg_Open_F(file_name, &
                 mode,      &
                 file_id,   &
                 error)

  if (error .ne. 0) then
    print *, "# Failed to open the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

  end subroutine
