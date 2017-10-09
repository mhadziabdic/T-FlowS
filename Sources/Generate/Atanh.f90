!======================================================================!
  real function Atanh(x)
!----------------------------------------------------------------------!
!   Calculates inverse of hyperbolic tangens.                          !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  real :: x
!======================================================================!

  if(x  > 1.0) then
    write(*,*) 'Error message from atanh: bad argument'
    stop
  end if 

  atanh=log( sqrt( (1.0+x)/(1.0-x) ) )

  end function Atanh