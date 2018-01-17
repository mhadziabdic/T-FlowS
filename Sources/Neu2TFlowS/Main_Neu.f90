!==============================================================================!
  program Neu2TFlowS 
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: grid     ! grid to be converted
  integer         :: c, n, s
!==============================================================================!

  call Logo_Neu

  write(*,*) '#======================================================'
  write(*,*) '# Enter the Fluent''s (*.NEU) file name (without ext.):'
  write(*,*) '#------------------------------------------------------'
  read(*,*) name

  call Load_Neu        (grid)
  call Grid_Topology   (grid)
  call Find_Faces      (grid)
  call Compute_Geometry(grid)
  call Connect_Domains (grid)

  do n=1,grid % n_nodes
    NewN(n) = n 
  end do  
  do c=-grid % n_bnd_cells,grid % n_cells
    NewC(c) = c 
  end do  
  do s=1,grid % n_faces 
    NewS(s) = s
  end do  

  call Save_Gmv_Cells(grid, 0,            &
                      grid % n_nodes,     &
                      grid % n_cells,     &
                      grid % n_faces,     &
                      grid % n_bnd_cells)

  call Save_Gmv_Faces(grid, 0,            &
                      grid % n_nodes)        ! save grid for checking b.c. 

  call Save_Shadows  (grid, 0,            &
                      grid % n_cells)        ! save shadows 

  call Save_Cns_Geo  (grid, 0,             &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0, 0) 

  ! Save links for checking
  call Save_Gmv_Links(grid, 0,             &
                      grid % n_nodes,      &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0)

  ! Create output for Fluent
  NewC(-grid % n_bnd_cells-1) = -grid % n_bnd_cells-1
  call Save_Cas(grid, 0,                        &
                grid % n_nodes,                 &
                grid % n_cells,                 &
                grid % n_faces + grid % n_sh,   &
                grid % n_bnd_cells)      ! save grid for postprocessing
                                         ! with Fluent

 ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)     ! save grid for postprocessing

  ! Create 1D file (used for channel or pipe flow) 
  call Probe_1D_Nodes(grid)

  ! Make .eps figures
  write(*,*) '# Making three .eps cuts through the domain.'
  call Save_Eps_Cut(grid, grid % dy, grid % dz, 'x')
  call Save_Eps_Cut(grid, grid % dz, grid % dx, 'y')
  call Save_Eps_Cut(grid, grid % dx, grid % dy, 'z')
 
  write(*,*) '# Making a 3D shaded .eps figure of the domain.'
  call Save_Eps_Whole(grid, grid % n_sh)   ! draw the domain with shadows

  end program
