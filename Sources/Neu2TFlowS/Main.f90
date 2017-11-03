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

  call Logo

  ! Test the precision
  open(90,FORM='unformatted',file='Neu2FlowS.real');
  write(90) 3.1451592
  close(90)  

  write(*,*) '#======================================================'
  write(*,*) '# Enter the Fluent''s (*.NEU) file name (without ext.):'
  write(*,*) '#------------------------------------------------------'
  read(*,*) name

  call Load_Neu(grid)
  call Grid_Topology(grid)
  call Find_Faces(grid)
  call Compute_Geometry(grid)
  
  call Connect_Domains(grid)

  do n=1,NN
    NewN(n) = n 
  end do  
  do c=-NbC,NC
    NewC(c) = c 
  end do  
  do s=1,NS 
    NewS(s) = s
  end do  

  call Save_Gmv_Cells(grid, 0, NN, NC, NS, NbC)
  call Save_Gmv_Faces(grid, 0, NN, NC)  ! save grid for checking b.c. 
  call Save_Shadows  (grid, 0, NN, NC)             ! save shadows 

  call Save_Cns_Geo(grid, 0, NC, NS, NBC, 0, 0) 

  ! Save links for checking
  call Save_Gmv_Links(grid, 0, NN, NC, NS, NbC, 0)

  ! Create output for Fluent
  NewC(-NBC-1) = -NBC-1
  call Save_Cas(grid, 0, NN, NC, NS+NSsh, NBC)  ! save grid for postprocessing
                                                ! with Fluent

  ! Create 1D file (used for channel or pipe flow) 
  call Probe_1D_Nodes(grid)

  ! Make .eps figures
  write(*,*) '# Making three .eps cuts through the domain.'
  call Save_Eps_Cut(grid, Dy, Dz, 'x')
  call Save_Eps_Cut(grid, Dz, Dx, 'y')
  call Save_Eps_Cut(grid, Dx, Dy, 'z')
 
  write(*,*) '# Making a 3D shaded .eps figure of the domain.'
  call Save_Eps_Whole(grid, NSsh)   ! Draw the domain with shadows

  end program
