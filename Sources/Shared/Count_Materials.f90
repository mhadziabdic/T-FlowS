!==============================================================================!
  subroutine Count_Materials()
!------------------------------------------------------------------------------!
!   Purpose: Counts all the materials in the grid.                             !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
!==============================================================================!

  Mater = .false.
  do c=1,NC  
    Mater(material(c)) = .true.
  end do

  Nmat = 0
  do c=1,1024
    if( Mater(c) ) Nmat = Nmat + 1
  end do

  write(*,*) 'Number of materials: ', Nmat
 
  end subroutine Count_Materials
