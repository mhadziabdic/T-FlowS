include '../Shared/To_Lower_Case.f90'
include '../Shared/To_Upper_Case.f90'
include 'Tokenizer.f90'

!==============================================================================!
  program Command_To_Control
!------------------------------------------------------------------------------!
!   Reads command file and turns it into control file.                         !
!------------------------------------------------------------------------------!
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: problem_name
  character(len=80) :: name_in
  character(len=80) :: answer
  integer           :: it
!==============================================================================!

  !---------------------------!
  !   Open the command file   !
  !---------------------------!
  open(CMN_FILE, file='T-FlowS.cmn')    
  cmn_line_count = 0

  !-----------------------!
  !   Read problem name   !
  !-----------------------!
  print *, '# Input problem name:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)  
  read(line % tokens(1), '(A80)')  problem_name

  ! Print entry for the control file
  if(answer .ne. 'SKIP') then
    print *, 'PROBLEM_NAME ', problem_name
  end if

  !-------------------------------!
  !   Reads type of the problem   !
  !-------------------------------!
  ! call Read_Problem           ! bad practice, should be avoided
  print *, '# Type of problem: '
  print *, '# CHANNEL          -> Channel flow'
  print *, '# PIPE             -> Pipe flow'
  print *, '# JET              -> Impinging jet flow'
  print *, '# TEST             -> Test Laplacian equation'
  print *, '# OTHER            -> All the other problems'
  print *, '# ROUGH            -> Problems with roughness'
  print *, '# HOT              -> Problems with temperature'
  print *, '# XHOM, YHOM, ZHOM -> Homogeneous directions'
  print *, '# TGV -> Taylor-Green Vortex test case'
  print *, '# BUOY -> Buoyancy flows (Automatically turns HOT on)' 
  print *, '# RB_CONV          -> Rayleigh-Barnard convection'
  print *, '# URANS            -> Unsteady RANS'
  print *, '# BACKSTEP         -> Backstep flow'
  do it = 1, line % n_tokens
    read(line % tokens(it),'(A8)')  answer
    call To_Lower_Case(answer)

    ! Print entry for the control file
    print *, 'PROBLEM_TYPE ', answer

    if(answer == 'rot') then
      print *, '# Angular velocity vector: '
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      ! Print entry for the control file
      print *, 'ANGULAR_VELOCITY_VECTOR ',  trim(line % tokens(1)),  &
                                            trim(line % tokens(2)),  &
                                            trim(line % tokens(3))

    else if(answer == 'BUOY') then
      print *, '# Gravitational constant in x, y and z directions: '
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      ! Print entry for the control file
      print *, 'GRAVITATIONAL_VECTOR ',  trim(line % tokens(1)),  &
                                         trim(line % tokens(2)),  &
                                         trim(line % tokens(3))
      print *, 'REFERENCE_TEMPERATURE ', trim(line % tokens(4))

    else if(answer == 'rough') then
      print *, '# Reading roughness coefficient Zo'
      call Tokenizer_Mod_Read_Line(CMN_FILE)

      ! Print entry for the control file
      print *, 'ROUGHNESS_COEFFICIENT ', trim(line % tokens(1))

    else
      print *, '# Error in input ! Exiting'
      stop
    endif
  end do

  !------------------------!
  !   Reads restart file   !
  !------------------------!
  ! call Load_Restart(grid, restar)
  print *, '# Input restart file name [skip cancels]:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1), '(A80)')  name_in
  answer=name_in
  call To_Upper_Case(answer) 

  ! Print entry for the control file
  if(answer .ne. 'SKIP') then
    print *, 'RESTART_FILE_NAME ', name_in
  end if

  !-----------------------!
  !   Reads T-Flows.cmn   !
  !-----------------------!
  ! call ReaCom(grid, restar)

  !--------------------------------------!
  !   Interpolate between diff. meshes   !
  !--------------------------------------!
  ! call Load_Ini(grid)

  !--------------------------------------------!
  !   Loading data from previous computation   !
  !--------------------------------------------!
  ! call Load_Restart_Ini(grid)

  !----------------------------!
  !   Close the command file   !
  !----------------------------!
  close(CMN_FILE)                

  end program
