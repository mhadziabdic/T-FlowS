!======================================================================!
  subroutine UserForce(var) 
!----------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: var ! 1 -> U,  2 -> V,  3 -> W 
!-------------------------------[Locals]-------------------------------!
  integer :: c
  real    :: buoyTemp, grav
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds bouyancy terms to the right hand side of velocity equation.   !
!                                                                      !
! Note:                                                                !
! ~~~~~                                                                !
!   Relies on two assumtions:                                          !
!   1. gravitational constant is assumed to be 1                       !
!   2. refference temperature is assumed to be 0                       !
!======================================================================!

  if (BUOY == NO) return ! XXXXX 5 Jul 2014

  if (var == 1) then
    grav = grav_x
  else if (var == 2) then
    grav = grav_y
  else if (var == 3) then
    grav = grav_z
  else
    grav = 0.0
  end if

  do c=1,NC
!    buoyTemp = 1.0/(273.15+T%n(c)) * grav*(T%n(c) - Tref)*volume(c)   
    buoyTemp = grav*(T%n(c) - Tref)*volume(c)   
    b(c) = b(c) - buoyTemp
  end do

  end subroutine UserForce
