!======================================================================!
  LOGICAL FUNCTION Approx(A,B,tol)
!----------------------------------------------------------------------!
!   Returns true if A~B, false otherwise.                              !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL          :: A,B 
  REAL,OPTIONAL :: tol
!-------------------------------[Locals]-------------------------------!
  REAL :: tolerance
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

  END FUNCTION Approx
