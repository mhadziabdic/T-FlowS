!==============================================================================!
  subroutine Load_Cns(grid)
!------------------------------------------------------------------------------!
!   Reads name.cns and first part of T-FlowS.cmn                               !
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: CMN_FILE
  use all_mod
  use pro_mod
  use par_mod,  only: this_proc
  use rans_mod, only: grav_x, grav_y, grav_z, Zo
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s, n, it
  character(len=80) :: name_in
  character(len=8)  :: answer
!==============================================================================!

  if(this_proc < 2) write(*,*) '# Input problem name:'
  call ReadC(CMN_FILE,inp,tn,ts,te)  
  read(inp(ts(1):te(1)), '(A80)')  name

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  call Name_File(this_proc, name_in, '.cns', len_trim('.cns'))  

  open(9, file=name_in,form='unformatted')
  if(this_proc < 2) write(*,*) '# Now reading the file:', name_in

  ! Number of cells, boundary cells and sides
  read(9) grid % n_cells                                 
  read(9) grid % n_boundary_cells
  read(9) grid % n_faces
  read(9) c   ! NSsh is not used in Processor. 

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Cells(grid, grid % n_boundary_cells, grid % n_cells) 
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces) 

  ! Number of materials and boundary conditions
  read(9) grid % n_materials
  read(9) grid % n_boundary_conditions

  allocate(grid % materials(grid % n_materials))
  allocate(grid % boundary_conditions(grid % n_boundary_conditions))

  ! Materials' and boundary conditions' keys
  do n=1,grid % n_materials
    read(9) grid % materials(n) % name
  end do
  do n=1,grid % n_boundary_conditions
    read(9) grid % boundary_conditions(n) % name
  end do

  ! Cell materials
  allocate (material(-grid % n_boundary_cells:grid % n_cells))
  read(9) (material(c), c =  1, grid % n_cells)        
  read(9) (material(c), c = -1,-grid % n_boundary_cells, -1) 

  ! Faces
  allocate (SideC(0:2,grid % n_faces))
  read(9) (SideC(0,s), s = 1, grid % n_faces)
  read(9) (SideC(1,s), s = 1, grid % n_faces)
  read(9) (SideC(2,s), s = 1, grid % n_faces)

  ! Boundary cells
  allocate (TypeBC(-grid % n_boundary_cells:grid % n_cells)); TypeBC=0 
  allocate (bcmark(-grid % n_boundary_cells:-1))
  read(9) (bcmark(c), c = -1, -grid % n_boundary_cells, -1) 

  ! Boundary copy cells
  allocate (CopyC(-grid % n_boundary_cells:-1))
  read(9) (CopyC(c), c = -1, -grid % n_boundary_cells, -1)   

  close(9)

  ! Type of the problem
  if(this_proc  < 2) then
    write(*,*) '# Type of problem: '
    write(*,*) '# CHANNEL          -> Channel flow'
    write(*,*) '# PIPE             -> Pipe flow'
    write(*,*) '# JET              -> Impinging jet flow'
    write(*,*) '# TEST             -> Test Laplacian equation'
    write(*,*) '# OTHER            -> All the other problems'
    write(*,*) '# ROUGH            -> Problems with roughness'
    write(*,*) '# HOT              -> Problems with temperature'
    write(*,*) '# XHOM, YHOM, ZHOM -> Homogeneous directions'
    write(*,*) '# TGV -> Taylor-Green Vortex test case'
    write(*,*) '# BUOY -> Buoyancy flows (Automatically turns HOT on)' 
  endif
  call ReadC(CMN_FILE,inp,tn,ts,te)
  do it=1,tn
    read(inp(ts(it):te(it)),'(A8)')  answer
    call To_Upper_Case(answer)
    if(answer == 'CHANNEL') then
      CHANNEL = YES
    else if(answer == 'PIPE') then
      PIPE = YES
    else if(answer == 'JET') then
      JET = YES
    else if(answer == 'TEST') then
      TEST = YES
    else if(answer == 'OTHER') then
      OTHER = YES
    else if(answer == 'HOT') then
      HOT = YES
    else if(answer == 'XHOM') then
      XHOM = YES
    else if(answer == 'YHOM') then
      YHOM = YES
    else if(answer == 'ZHOM') then
      ZHOM = YES
    else if(answer == 'TGV') then
      TGV = YES
    else if(answer == 'ROT') then
      ROT = YES
    else if(answer == 'BUOY') then
      BUOY = YES
    else if(answer == 'RB_CONV') then
      RB_CONV = YES
      BUOY = YES
      HOT  = YES
    else if(answer == 'BUDG') then
      BUDG = YES
    else if(answer == 'BACKSTEP') then
      BACKSTEP = YES
    else if(answer == 'URANS') then
      URANS = YES
    else if(answer == 'ROUGH') then
      ROUGH = YES
    else
      write(*,*) 'Error in input ! Exiting'
      stop
    endif
  end do

  if(ROUGH == YES) then
    if(this_proc < 2) write(*,*) '# Reading roughness coefficient Zo'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*) Zo
  endif

  ! Angular velocity vector
  if(ROT == YES) then
    if(this_proc  < 2)  &
    write(*,*) '# Angular velocity vector: '
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)  omegaX, omegaY, omegaZ
  end if

  ! Gravity
  if(BUOY == YES) then
    if(this_proc  < 2)  &
    write(*,*) '# Gravitational constant in x, y and z directions: '
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*) grav_x, grav_y, grav_z, Tref
  end if

  end subroutine Load_Cns
