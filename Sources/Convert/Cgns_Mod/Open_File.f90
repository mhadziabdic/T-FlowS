!==============================================================================!
  subroutine Cgns_Mod_Open_File
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer*8 :: error
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

  !  set initial values
  n_nodes = 0
  n_cells = 0
  n_hex   = 0
  n_pyr   = 0
  n_wed   = 0
  n_tet   = 0
  n_tri   = 0
  n_qua   = 0

  last_x = 0
  last_y = 0
  last_z = 0

  last_hex = 0
  last_pyr = 0
  last_wed = 0
  last_tet = 0
  last_tri  = 0
  last_qua = 0

  end subroutine
