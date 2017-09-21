!======================================================================!
  REAL FUNCTION Atanh(x)
!----------------------------------------------------------------------!
!   Calculates inverse of hyperbolic tangens.                          !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: x
!--------------------------------[CVS]---------------------------------!
!  $Id: Atanh.f90,v 1.1 2014/11/24 11:31:30 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Generate/Atanh.f90,v $  
!======================================================================!
  if(x  > 1.0) then
    write(*,*) 'Error message from atanh: bad argument'
    stop
  end if 

  atanh=log( sqrt( (1.0+x)/(1.0-x) ) )

  END FUNCTION Atanh
