!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Domain_Mod  ! domain as defined in ".dom" file.
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type) :: dom    ! domain to be used
  type(Grid_Type)   :: grid   ! grid which will be generated
  integer           :: c, s, n
!==============================================================================!

  ! Open with a logo
  call Logo_Gen

  call Load_Domain             (dom, grid)
  call Compute_Node_Coordinates(dom, grid)
  call Distribute_Regions      (dom, grid)
  call Connect_Blocks          (dom, grid)
  call Connect_Periodicity     (dom, grid)
  call Connect_Copy            (dom)

  ! From this point on, domain is not used anymore
  call Determine_Grid_Connectivity(grid, .false.)  ! trial run 
  call Compute_Grid_Geometry      (grid, .false.)
  call Smooth_Grid                (grid)
  call Refine_Grid                (grid)
  call Determine_Grid_Connectivity(grid, .true.) ! real run
  call Compute_Grid_Geometry      (grid, .true.)

  ! Prepare for saving
  do n = 1,grid % n_nodes
    NewN(n)=n
  end do
  do c = -grid % n_bnd_cells,grid % n_cells
    NewC(c)=c
  end do
  do s = 1,grid % n_faces
    NewS(s)=s
  end do

  ! Save the grid
  call Save_Gmv_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)     ! save grid for postprocessing

  call Save_Gmv_Faces(grid, 0,         &
                      grid % n_nodes)     ! save grid for checking b.c. 

  call Save_Shadows  (grid, 0,         &
                      grid % n_cells)     ! save shadows 

  call Save_Cns_Geo(grid, 0,                  &
                    grid % n_cells,           &
                    grid % n_faces,           &
                    grid % n_bnd_cells,  &
                    0, 0)  ! saved data for processing

  ! Save links for checking
  call Save_Gmv_Links(grid, 0,                  &
                      grid % n_nodes,           &
                      grid % n_cells,           &
                      grid % n_faces,           &
                      grid % n_bnd_cells,  &
                      0)

  ! Save the 1D probe (good for the channel flow)
  call Probe_1D_Nodes_Gen(grid)

  ! Save the 2D probe (good for the channel flow)
  call Probe_2D(grid)

  ! Create output for Fluent
  NewC(-grid % n_bnd_cells-1) = -grid % n_bnd_cells-1
  call Save_Cas(grid, 0,              &
                grid % n_nodes,       &
                grid % n_cells,       &
                grid % n_faces + grid % n_sh)  ! save grid for Fluent

  ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)     ! save grid for postprocessing

  ! Make eps figures
  call Save_Eps_Cut(grid, grid % dy, grid % dz, 'x') 
  call Save_Eps_Cut(grid, grid % dz, grid % dx, 'y') 
  call Save_Eps_Cut(grid, grid % dx, grid % dy, 'z') 

  call Save_Eps_Whole(grid, grid % n_sh)  ! draw the domain with shadows

  ! Write something on the screen
  call Print_Grid_Statistics(grid)

  end program
