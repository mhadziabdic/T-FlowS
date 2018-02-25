!==============================================================================!
  subroutine Read_Problem
!------------------------------------------------------------------------------!
!   Reads first part of T-FlowS.cmn. (BN: I can't see a good reason to keep    !
!   this separate from reading the rest of command file.)                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod, only: this_proc
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
    print *, '# Type of problem: '
    print *, '# CHANNEL          -> Channel flow'
    print *, '# PIPE             -> Pipe flow'
    print *, '# JET              -> Impinging jet flow'
    print *, '# TEST             -> Test Laplacian equation'
    print *, '# OTHER            -> All the other problems'
    print *, '# ROUGH            -> Problems with roughness'
    print *, '# XHOM, YHOM, ZHOM -> Homogeneous directions'
    print *, '# TGV -> Taylor-Green Vortex test case'
    print *, '# BUOY -> Buoyancy flows' 
  endif
  ! call Tokenizer_Mod_Read_Line(CMN_FILE)
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
!     BUOY = YES
    else if(answer == 'RB_CONV') then
      RB_CONV = YES
!     BUOY = YES
    else if(answer == 'BUDG') then
      BUDG = YES
    else if(answer == 'BACKSTEP') then
      BACKSTEP = YES
    else if(answer == 'ROUGH') then
      ROUGH = YES
    else
      print *, 'Error in input ! Exiting'
      stop
    endif
  end do

  if(ROUGH == YES) then
    if(this_proc < 2) print *, '# Reading roughness coefficient Zo'
!   call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), *) Zo
  endif

  ! Angular velocity vector
  if(ROT == YES) then
    if(this_proc  < 2)  &
    print *, '# Angular velocity vector: '
!   call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), *)  omega_x
    read(line % tokens(2), *)  omega_y
    read(line % tokens(3), *)  omega_z
  end if

! ! Gravity
! if(BUOY == YES) then
!   if(this_proc  < 2)  &
!   print *, '# Gravitational constant in x, y and z directions: '
!   call Tokenizer_Mod_Read_Line(CMN_FILE)
!   read(line % tokens(1), *) grav_x
!   read(line % tokens(2), *) grav_y
!   read(line % tokens(3), *) grav_z
!   read(line % tokens(4), *) Tref
! end if

  end subroutine
