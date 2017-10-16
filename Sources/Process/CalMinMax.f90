!======================================================================!
  subroutine CalMinMax(PHI)
!----------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  real    :: PHI(-NbC:NC) 
!-------------------------------[Locals]-------------------------------!
  integer :: c1, c2, s
!======================================================================!

  PHImax = PHI 
  PHImin = PHI 

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if( (c2>0) .or. (c2<0 .and. TypeBC(c2)==BUFFER) ) then
      PHImax(c1) = max(PHImax(c1),PHI(c2))
      PHImin(c1) = min(PHImin(c1),PHI(c2))
      PHImax(c2) = max(PHImax(c2),PHI(c1))
      PHImin(c2) = min(PHImin(c2),PHI(c1))
    end if

  end do

  call exchng(PHImax)
  call exchng(PHImin)

  RETURN 

  end subroutine CalMinMax
