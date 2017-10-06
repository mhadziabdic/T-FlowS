!======================================================================!
  subroutine CorBad(PHIi)
!----------------------------------------------------------------------!
!   Corrects the pressure gradients in the cells where they cannot     !
!   be computed, the so called "bad" cells.                            !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  real :: PHIi(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  integer :: c, c1, c2, s
!======================================================================!

  do c=1,NC
    if(BadForG(c)) then
      PHIi(c) = 0.0
    end if
  end do 

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
     
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then
      if(BadForG(c1)) PHIi(c1) = PHIi(c1) + 0.5*PHIi(c2) 
      if(BadForG(c2)) PHIi(c2) = PHIi(c2) + 0.5*PHIi(c1) 
    end if
  end do

  call Exchng(PHIi)

  end subroutine CorBad 
