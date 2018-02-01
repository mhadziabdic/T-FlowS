!==============================================================================!
  subroutine Cgns_Mod_Print_Base_Info
!------------------------------------------------------------------------------!
!   Reads main info from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: base_name
  integer           :: phys_dim
  integer           :: cell_dim
!==============================================================================!

  ! Read CGNS base information
  call Cg_Base_Read_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      base_name, & ! name of the base
                      cell_dim,         & ! cell dimensions (3->volume cell)
                      phys_dim,         & ! number of coordinates to create vector
                      ier)                ! error status

    if (ier .ne. 0) then
      print *, "# Failed to get base info"
      call Cg_Error_Exit_F()
    endif

  print *, "# Base Name: ",      base_name
  print *, "# Cell dimension: ", cell_dim
  print *, "# Phys dimension: ", phys_dim

  end subroutine
