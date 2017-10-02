!======================================================================!
  SUBROUTINE CouMat()
!----------------------------------------------------------------------!
! Purpose: Counts all the materials in the grid.                       !
! ~~~~~~~~                                                             !
!------------------------------[Modules]-------------------------------!
  USE all_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c 
!======================================================================!

  Mater = .FALSE.
  do c=1,NC  
    Mater(material(c)) = .TRUE.
  end do

  Nmat = 0
  do c=1,1024
    if( Mater(c) ) Nmat = Nmat + 1
  end do

  write(*,*) 'Number of materials: ', Nmat
 
  END SUBROUTINE CouMat
