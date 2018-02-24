!==============================================================================!
  program Neu2TFlowS 
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type)   :: grid     ! grid to be converted
  integer           :: c, n, s, l
  character(len=80) :: file_name, file_name_up, extension
!==============================================================================!

  call Logo_Con

  print *, '#==============================================================='
  print *, '# Enter the name of the mesh file you are importing (with ext.):'
  print *, '#---------------------------------------------------------------'
  read(*,*) file_name

  !-----------------------------------------------------!
  !   Analyze the extension to figure the file format   !
  !-----------------------------------------------------!
  file_name_up = file_name
  call To_Upper_Case(file_name_up)
 
  l = len_trim(file_name)
  if( file_name_up(l-2:l) == 'NEU' ) then
    print *, '# Based on the extension, you are' // &
             ' reading Gambit''s neutral file format'
    problem_name = file_name(1:l-4)
    extension = file_name_up(l-2:l)
  else if( file_name_up(l-3:l) == 'CGNS' ) then
    print *, '# Based on the extension, you are' // &
             ' reading CGNS file format'
    problem_name = file_name(1:l-5)
    extension = file_name_up(l-3:l)
  else
    print *, '# Unrecognized input file format; exiting!'
    stop
  end if

  !----------------------------------------!
  !   Read the file and start conversion   !
  !----------------------------------------!
  if (extension == 'NEU') then
    call Load_Neu (grid)
  end if
  if (extension == 'CGNS') then
    call Load_Cgns(grid)
  end if

  call Grid_Topology     (grid)
  call Find_Faces        (grid)
  call Calculate_Geometry(grid)
  call Connect_Domains (grid)

  ! Prepare for saving
  do n=1,grid % n_nodes
    NewN(n) = n 
  end do  
  do c=-grid % n_bnd_cells,grid % n_cells
    NewC(c) = c 
  end do  
  do s=1,grid % n_faces 
    NewS(s) = s
  end do  

  !-------------------------------!
  !   Save files for processing   !
  !-------------------------------!
  call Save_Cns_Geo  (grid, 0,             &
                      grid % n_nodes,      &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0, 0) 

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

 ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)

  ! Save links for checking
  call Save_Vtu_Links(grid, 0,             &
                      grid % n_nodes,      &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0)

  ! I believe these are needed for Fluent only
  call Save_Shadows(grid, 0,            &
                    grid % n_cells)

  ! Create 1D file (used for channel or pipe flow) 
  call Probe_1D_Nodes(grid)

  end program
