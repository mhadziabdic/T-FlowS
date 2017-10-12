!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                    ________  ______                                  !
!                   |        ||      \                                 !
!                   `--.  .--'|  ,-.  \________  ___                   !
!                      |  |___|  |_/  /  _  \  \/  /                   !
!                      |  |___|      /  _____\    /                    !
!                      |  |   |  |\  \  \____/    \                    !
!                      |__|   |__| \__\_____/__/\__\                   !
!                                                                      !
!                                                                      !
!           BLOCK-STRUCTURED 3D HEXAHEDRAL MESH GENERATOR              !
!                                +                                     !
!                  UNSTRUCTURED CELL REFINEMENT                        !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!                                     Bojan Niceno                     !
!                                     Delft University of Technology   !
!                                     Section Heat Transfer            !
!                                     niceno@duttwta.wt.tn.tudelft.nl  !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!  Block structure:                                                    !
!                                                                      !
!                               b6                                     !
!                               |                                      !
!                         8-----|---------6                            !
!                        /|     |        /|                            !
!                       / |     +   b3  / |                            !
!                      /  |        /   /  |                            !
!                     7---------------5   |                            !
!                     |   |     |/    |   |                            !
!              b4---- | +-------o-----| +-------b2                     !
!                     |   |    /|     |   |                            !
!                     |   4---/-|-----|---2                            !
!                     |  /   /  |     |  /                             !
!                     | /   b5  +     | /                              !
!                     |/              |/                               !
!                     3---------------1                                !
!                               |                                      !
!                               b1                                     !
!  Faces are defined as:
!
!     I   :
!     II  :
!     III :
!     IV  :
!     V   :
!     VI  :
!
!  Local coordinate directions:                                        !
!                                                                      !
!     i: 1 -> 2                                                        !
!     j: 1 -> 3                                                        !
!     k: 1 -> 5                                                        !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!    - can't handle domains with less then 3x3x3 cells properly        !
!    - nodes of a cell (and block) are deonoted with numbers 1 - 8     !
!    - neighbouring cells are denoted with c1 - c6                     !
!    - local coordinare directions (for blocks) are defined with:      !
!                                                                      !
!======================================================================!
  PROGRAM Generator
!----------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement. *
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c,s,n
!======================================================================!

  ! Test the precision
  open(90,FORM='UNFORMATTED',FILE='Generator.real'); 
  write(90) 3.1451592
  close(90)

  call Logo

  call Load_Domain
  call Compute_Node_Coordinates
  call Connect_Blocks
  call PeriBC
  call CopyBC

  call TopSys(.false.)  ! trial run 
  call Compute_Grid_Geometry(.false.)
  call Smooth 

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

  ! Count the materials in the grid
  call Count_Materials

  ! Save the grid
  call Save_Gmv_Mesh(0, NN, NC)            ! save grid for postprocessing
  call GeoSav(0, NC, NS, NBC, 0, 0) ! saved data for processing

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
  call Save_Eps_Cut(y_node,z_node,x_node,Dy,Dz,'x') 
  call Save_Eps_Cut(z_node,x_node,y_node,Dz,Dx,'y') 
  call Save_Eps_Cut(x_node,y_node,z_node,Dx,Dy,'z') 

  call Save_Eps_Whole(NSsh)  ! draw the domain with shadows

  ! Write something on the screen
  call Print_Grid_Statistics

  end PROGRAM Generator 
