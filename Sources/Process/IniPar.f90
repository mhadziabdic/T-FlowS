!======================================================================!
  subroutine IniPar 
!----------------------------------------------------------------------!
!   Initialize various solver parameters.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use Flow_Mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: i
!======================================================================!

  !---------------------------!
  !   Initialize parameters   !
  !---------------------------!

  ! for wall functions
  CHANNEL = NO     ! YES, NO
  BACKSTEP= NO     ! YES, NO
  RB_CONV = NO     ! YES, NO
  PIPE    = NO     ! YES, NO
  JET     = NO     ! YES, NO
  ROUGH   = NO     ! YES, NO
  TGV    = NO      ! YES, NO
  TEST   = NO      ! YES, NO
  OTHER  = NO      ! YES, NO
  SGDH   = NO      ! YES, NO
  GGDH   = NO      ! YES, NO
  ROT    = NO      ! YES, NO
  HOTini = NO      ! YES, NO
  SHAKE  = NO      ! YES, NO
  BUDG   = NO
  SHAKE_PER = -10
  SHAKE_INT = -10
  XHOM   = NO
  YHOM   = NO
  ZHOM   = NO
  PER_BC = YES

  ! Set the default monitoring point
  Nmon=1
  do i=1,Nmon
    Cm(i) = 1
  end do

  end subroutine
