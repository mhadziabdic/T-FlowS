!==============================================================================!
  PROGRAM Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: c,s,n
!==============================================================================!

  ! Test the precision
  open(90,FORM='unformatted',file='Generator.real'); 
  write(90) 3.1451592
  close(90)

  call Logo

  call Load_Domain
  call Compute_Node_Coordinates
  call Distribute_Regions
  call Connect_Blocks
  call Connect_Periodicity
  call CopyBC

  call TopSys(.false.)  ! trial run 
  call Compute_Grid_Geometry(.false.)
  call Smooth_Grid

  call Refine_Grid

  call TopSys(.true.) ! real run
  call Compute_Grid_Geometry(.true.)

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
  call Save_Gmv_Grid(0, NN, NC)            ! save grid for postprocessing
  call Save_Cns_Geo(0, NC, NS, NBC, 0, 0)  ! saved data for processing

  call Save_Gmv_Links(0, NN, NC, NS, NbC, 0)

  ! Save the 1D probe (good for the channel flow)
  call Probe_1D_Nodes_Gen

  ! Save the 2D probe (good for the channel flow)
  call Probe_2D

  ! Create output for Fluent
  NewC(-NBC-1) = -NBC-1
  call Save_Cas(0, NN, NC, NS+NSsh) ! save grid for postprocessing
                                    ! with Fluent
  ! Make eps figures
  call Save_Eps_Cut(Dy,Dz,'x') 
  call Save_Eps_Cut(Dz,Dx,'y') 
  call Save_Eps_Cut(Dx,Dy,'z') 

  call Save_Eps_Whole(NSsh)  ! draw the domain with shadows

  ! Write something on the screen
  call Print_Grid_Statistics

  end PROGRAM Generator 
