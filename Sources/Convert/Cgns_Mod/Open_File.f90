!==============================================================================!
  subroutine Cgns_Mod_Open_File
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  print *, "# Reading the file:", trim(file_name)

  ! Open a CGNS file
  call Cg_Open_F(file_name,    &  ! file name
                 CG_MODE_READ, &  ! open for read
                 file_id,      &  ! cgns file index number
                 ier)             ! error status

  if (ier .ne. 0) then
    print *, "# Failed to read the file: ", trim(file_name)
    call Cg_Error_Exit_F()
  endif

  !  set initial values
  n_nodes = 0
  n_cells = 0
  n_hexa  = 0
  n_pyra  = 0
  n_pris  = 0
  n_tetr  = 0
  n_tria  = 0
  n_quad  = 0

  last_x = 0
  last_y = 0
  last_z = 0

  last_hexa = 0
  last_pyra = 0
  last_pris = 0
  last_tetr = 0
  last_tria = 0
  last_quad = 0

  bc_id = 0

  end subroutine
