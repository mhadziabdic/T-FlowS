!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Domain_Mod  ! domain as defined in ".d" file.
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type) :: dom    ! domain to be used
  type(Grid_Type)   :: grid   ! grid which will be generated
  integer           :: c, s, n
!==============================================================================!

  ! Test the precision
  open(90,FORM='unformatted',file='Generator.real'); 
  write(90) 3.1451592
  close(90)

  ! Open with a logo
  call Logo

  call Load_Domain(dom, grid)
  call Compute_Node_Coordinates(dom, grid)
  call Distribute_Regions(dom, grid)
  call Connect_Blocks(dom, grid)
  call Connect_Periodicity(dom, grid)
  call Connect_Copy(dom)

  ! From this point on, domain is not used anymore
  call Determine_Grid_Connectivity(grid, .false.)  ! trial run 
  call Compute_Grid_Geometry(grid, .false.)
  call Smooth_Grid(grid)

  call Refine_Grid(grid)

  call Determine_Grid_Connectivity(grid, .true.) ! real run
  call Compute_Grid_Geometry(grid, .true.)

  ! Prepare for saving
  do n=1,NN
    NewN(n)=n
  end do
  do c=-NBC,NC
    NewC(c)=c
  end do
  do s=1,NS
    NewS(s)=s
  end do

  ! Save the grid
  call Save_Gmv_Grid(grid, 0, NN, NC)            ! save grid for postprocessing
  call Save_Cns_Geo(grid, 0, NC, NS, NBC, 0, 0)  ! saved data for processing

  ! Save links for checking
  call Save_Gmv_Links(grid, 0, NN, NC, NS, NbC, 0)

  ! Save the 1D probe (good for the channel flow)
  call Probe_1D_Nodes_Gen(grid)

  ! Save the 2D probe (good for the channel flow)
  call Probe_2D(grid)

  ! Create output for Fluent
  NewC(-NBC-1) = -NBC-1
  call Save_Cas(grid, 0, NN, NC, NS+NSsh) ! save grid for postprocessing
                                    ! with Fluent
  ! Make eps figures
  call Save_Eps_Cut(grid, Dy,Dz,'x') 
  call Save_Eps_Cut(grid, Dz,Dx,'y') 
  call Save_Eps_Cut(grid, Dx,Dy,'z') 

  call Save_Eps_Whole(grid, NSsh)  ! draw the domain with shadows

  ! Write something on the screen
  call Print_Grid_Statistics(grid)

  end program
