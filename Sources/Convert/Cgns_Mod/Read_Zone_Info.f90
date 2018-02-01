!==============================================================================!
  subroutine Cgns_Mod_Read_Zone_Info
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Read zone information
  call Cg_Zone_Read_F(file_id,   & ! cgns file index number
                      base_id,   & ! base index number
                      zone_id,   & ! zone index number
                      zone_name, & ! name of the zone
                      mesh_info, & ! n_nodes, n_cells, n_b_nodes(if sorted)
                      ier)         ! error status

  if (ier .ne. 0) then
    print *, "#     Failed read zone info"
    call Cg_Error_Exit_F()
  endif

  ! total nodes and cells
  n_nodes = n_nodes + mesh_info(1)
  n_cells = n_cells + mesh_info(2)

  print *, "#     Zone index: ",                zone_id
  print *, "#     Zone name: ",                 zone_name
  print *, "#     Nodes: ",                     mesh_info(1)
  print *, "#     Cells: ",                     mesh_info(2)
  print *, "#     Boundary nodes(if sorted): ", mesh_info(3)

  if (mesh_info(3) .ne. 0) then
    print *, "#     B.C. nodes != 0 -> Unsupported"
    stop
  endif

  end subroutine
