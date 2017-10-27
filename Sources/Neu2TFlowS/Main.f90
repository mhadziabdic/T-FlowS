!==============================================================================!
  program Neu2TFlowS 
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, s
!==============================================================================!

  call Logo

  ! Test the precision
  open(90,FORM='unformatted',file='Neu2FlowS.real');
  write(90) 3.1451592
  close(90)  

  write(*,*) '#======================================================'
  write(*,*) '# Enter the Fluent''s (*.NEU) file name (without ext.):'
  write(*,*) '#------------------------------------------------------'
  read(*,*) name

  call Load_Neu
  call Grid_Topology
  call Find_Sides
  call Compute_Geometry
  
  call Connect_Domains

  do n=1,NN
    NewN(n) = n 
  end do  
  do c=-NbC,NC
    NewC(c) = c 
  end do  
  do s=1,NS 
    NewS(s) = s
  end do  

  ! Count all materials
  call Count_Materials

  call Save_Gmv_Grid(0, NN, NC, NS, NbC)
  call Save_Cns_Geo(0, NC, NS, NBC, 0, 0) 
  call Save_Gmv_Links(0, NN, NC, NS, NbC, 0)

  ! Create output for Fluent
  NewC(-NBC-1) = -NBC-1
! call Save_Cas(0, NN, NC, NS+NSsh, NBC) ! save grid for postprocessing
                                         ! with Fluent

  ! Create 1D file (used for channel or pipe flow) 
  call Probe_1D_Nodes

  ! Make .eps figures
  write(*,*) '# Making three .eps cuts through the domain.'
  call Save_Eps_Cut(Dy, Dz, 'x')
  call Save_Eps_Cut(Dz, Dx, 'y')
  call Save_Eps_Cut(Dx, Dy, 'z')
 
  write(*,*) '# Making a 3D shaded .eps figure of the domain.'
  call Save_Eps_Whole(NSsh)   ! Draw the domain with shadows

  end program Neu2TFlowS 
