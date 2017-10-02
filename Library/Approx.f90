!======================================================================!
  logical function Approx(A,B,tol)
!----------------------------------------------------------------------!
!   Returns true if A~B, false otherwise.                              !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  real          :: A,B 
  real,OPTIONAL :: tol
!-------------------------------[Locals]-------------------------------!
  real :: tolerance
!======================================================================!

  if( .not. present(tol) ) then
    tolerance = 1.e-9
  else
    tolerance = tol
  end if

  if( (A  < (B + tolerance)) .and. (A  > (B - tolerance)) ) then
    approx = .TRUE.
  else
    approx = .FALSE.
  end if

  end function Approx
