!======================================================================!
  subroutine UserZero(var) 
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds bouyancy terms to the right hand side of velocity equation    !
!   for zero Pr cavity L x H = 4 x 1.                                  !
!                                                                      ! 
! Note:                                                                !
! ~~~~~                                                                !
!   It should be used with zero.* problem in Test directory.           !  
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: var ! 1 -> U,  2 -> V,  3 -> W 
!-------------------------------[Locals]-------------------------------!
  integer :: c
!======================================================================!

  if(var == 3) then  ! only for W velocity component
    do c=1,NC
      b(c)=b(c) + 0.25*(xc(c)-2.0)*volume(c)
    end do
  end if

  end subroutine UserZero
