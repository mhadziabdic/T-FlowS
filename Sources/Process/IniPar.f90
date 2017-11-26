!======================================================================!
  subroutine IniPar 
!----------------------------------------------------------------------!
!   Initialize various solver parameters.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Parameters_Mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: i
!======================================================================!

  !---------------------------!
  !   Initialize parameters   !
  !---------------------------!

  ! for wall functions
  MODE       = 28

  ALGOR  = FRACT   ! FRACT, SIMPLE
  INERT  = LIN     ! LIN, PAR 
  CONVEC = AB      ! AB, CN, FI
  CROSS  = AB      ! AB, CN, FI
  DIFFUS = CN      ! AB, CN, FI
  SIMULA = DNS     ! LES, DNS  
  CHANNEL = NO     ! YES, NO
  BACKSTEP= NO     ! YES, NO
  RB_CONV = NO     ! YES, NO
  PIPE    = NO     ! YES, NO
  JET     = NO     ! YES, NO
  ROUGH   = NO     ! YES, NO
  TGV    = NO      ! YES, NO
  TEST   = NO      ! YES, NO
  OTHER  = NO      ! YES, NO
  HOT    = NO      ! YES, NO
  SGDH   = NO      ! YES, NO
  GGDH   = NO      ! YES, NO
  ROT    = NO      ! YES, NO
  HOTini = NO      ! YES, NO
  SHAKE  = NO      ! YES, NO
  BLEND  = NO      ! YES, NO, QUICK,, LUDS, MINMOD, SMART, AVL_SMART
  BUOY   = NO
  URANS  = NO
  BUDG   = NO
  SHAKE_PER = -10
  SHAKE_INT = -10
  XHOM   = NO
  YHOM   = NO
  ZHOM   = NO
  PER_BC = YES

  ! Under relaxation factors
  U % URF      = 1.0
  P % URF      = 1.0
  T % URF      = 1.0
  Kin % URF    = 1.0
  Eps % URF    = 1.0
  v_2 % URF    = 1.0
  VIS % URF    = 1.0
  URFC         = 0.9
  
  Ndyn         = 0
 
  ! Algorythm and solver tolerances
  SIMTol       = 1.e-4
  U % STol     = 1.e-6 
  T % STol     = 1.e-6 
  PP % STol    = 1.e-6
  Kin % STol   = 1.e-6
  Eps % STol   = 1.e-6
  v_2 % STol   = 1.e-6
  f22 % STol   = 1.e-6
  VIS % STol   = 1.e-6

  ! Set the default monitoring point
  Nmon=1
  do i=1,Nmon
    Cm(i) = 1
  end do

  end subroutine
