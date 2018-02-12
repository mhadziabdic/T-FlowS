!==============================================================================!
  subroutine Cgns_Mod_Open_File
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  print *, "# Reading the file:", trim(file_name)

  ! Open a CGNS file
  call Cg_Open_F(file_name,    &
                 CG_MODE_READ, &
                 file_id,      &
                 error)

  if (error .ne. 0) then
    print *, "# Failed to read the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

  end subroutine
