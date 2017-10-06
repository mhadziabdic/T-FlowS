!======================================================================!
  subroutine CnsLoa
!----------------------------------------------------------------------!
! Reads NAME.cns and first part of T-FlowS.cmn                         !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  use allp_mod, only: CMN_FILE
  use all_mod
  use pro_mod
  use par_mod,  only: this
  use rans_mod, only: grav_x, grav_y, grav_z, Zo
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer      :: c, s, it
  character*80 :: nameIn
  character*8  :: answer
!======================================================================!

  if(this < 2) write(*,*) '# Input problem name:'
  call ReadC(CMN_FILE,inp,tn,ts,te)  
  read(inp(ts(1):te(1)), '(A80)')  name

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!       Read the file with the      !
!     connections between cells     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call NamFil(THIS, nameIn, '.cns', len_trim('.cns'))  

  open(9, FILE=nameIn,FORM='UNFORMATTED')
  if(this < 2) write(*,*) '# Now reading the file:', nameIn

!///// number of cells, boundary cells and sides
  read(9) NC                                 
  read(9) NbC
  read(9) NS
  read(9) c   ! NSsh is not used in Processor. 
  read(9) Nmat

!///// cell materials
  allocate (material(-NbC:NC))
  read(9) (material(c), c=1,NC)        
  read(9) (material(c), c=-1,-NBC,-1) 

!///// sides
  allocate (SideC(0:2,NS))
  read(9) (SideC(0,s), s=1,NS)
  read(9) (SideC(1,s), s=1,NS)
  read(9) (SideC(2,s), s=1,NS)

!///// boundary cells
  allocate (TypeBC(-NbC:NC)); TypeBC=0 
  allocate (bcmark(-NbC:-1))
  read(9) (bcmark(c), c=-1,-NbC, -1) 

!///// boundary copy cells
  allocate (CopyC(-NbC:-1))
  read(9) (CopyC(c), c=-1,-NbC, -1)   

  close(9)

!----- Type of the problem
  if(this  < 2) then
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
    call ToUppr(answer)
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
    if(this < 2) write(*,*) '# Reading roughness coefficient Zo'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*) Zo
  endif

!----- Angular velocity vector
  if(ROT == YES) then
    if(this  < 2)  &
    write(*,*) '# Angular velocity vector: '
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)  omegaX, omegaY, omegaZ
  end if

!----- Gravity
  if(BUOY == YES) then
    if(this  < 2)  &
    write(*,*) '# Gravitational constant in x, y and z directions: '
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*) grav_x, grav_y, grav_z, Tref
  end if

  end subroutine CnsLoa
