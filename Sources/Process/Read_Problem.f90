!==============================================================================!
  subroutine Read_Problem
!------------------------------------------------------------------------------!
!   Reads first part of T-FlowS.cmn. (BN: I can't see a good reason to keep    !
!   this separate from reading the rest of command file.)                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: CMN_FILE
  use all_mod
  use pro_mod
  use par_mod,  only: this_proc
  use rans_mod, only: grav_x, grav_y, grav_z, Zo
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: it
  character(len=8)  :: answer
!==============================================================================!

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
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  do it = 1, line % n_tokens
    read(line % tokens(it),'(A8)')  answer
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
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), *) Zo
  endif

  ! Angular velocity vector
  if(ROT == YES) then
    if(this_proc  < 2)  &
    write(*,*) '# Angular velocity vector: '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), *)  omegaX
    read(line % tokens(2), *)  omegaY
    read(line % tokens(3), *)  omegaZ
  end if

  ! Gravity
  if(BUOY == YES) then
    if(this_proc  < 2)  &
    write(*,*) '# Gravitational constant in x, y and z directions: '
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), *) grav_x
    read(line % tokens(2), *) grav_y
    read(line % tokens(3), *) grav_z
    read(line % tokens(4), *) Tref
  end if

  end subroutine
