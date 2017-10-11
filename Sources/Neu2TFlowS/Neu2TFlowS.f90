!======================================================================!
  program Neu2TFlowS 
!------------------------------[Modules]-------------------------------!
  use all_mod 
  use gen_mod 
!----------------------------------------------------------------------! 
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, n, s
  integer :: i, typ, cnt
!======================================================================!

  call Logo

!---- Test the precision
  open(90,FORM='UNFORMATTED',FILE='Neu2FlowS.real');
  write(90) 3.1451592
  close(90)  

  call ReadFluentNeu
  call TopolM
  call FindSides
  call Calc4
  
  call Connect3

  do n=1,NN
    NewN(n) = n 
  end do  
  do c=-NbC,NC
    NewC(c) = c 
  end do  
  do s=1,NS 
    NewS(s) = s
  end do  

!---- caount all the materials
  call CouMat

  call GenSav(0, NN, NC, NS, NbC)
  call GeoSav(0, NC, NS, NBC, 0, 0) 
  call TestLn(0, NN, NC, NS, NbC, 0)

!---- create output for Fluent
  NewC(-NBC-1) = -NBC-1
! call CasSav(0, NN, NC, NS+NSsh, NBC) ! Save grid for postprocessing
                                       ! with Fluent

!---- create 1D file (used for channel or pipe flow) 
  call Probe_1D_Nodes

!---- make eps figures
  write(*,*) 'Making three .eps cuts through the domain.'
  call EpsSav(y_node, z_node, x_node, Dy, Dz, 'x')
  call EpsSav(z_node, x_node, y_node, Dz, Dx, 'y')
  call EpsSav(x_node, y_node, z_node, Dx, Dy, 'z')
 
  write(*,*) 'Making a 3D shaded .eps figure of the domain.'
  call EpsWho(NSsh)  ! Draw the domain with shadows

  end program Neu2TFlowS 
