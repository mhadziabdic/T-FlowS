!======================================================================!
  LOGICAL FUNCTION IsTwin(n1, n2) 
!----------------------------------------------------------------------!
!   Checks if the nodes are twins, i.e. are they shared on periodicity !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: n1, n2
!-------------------------------[Locals]-------------------------------!
  INTEGER :: n, check
!--------------------------------[CVS]---------------------------------!
!  $Id: IsTwin.f90,v 1.1 2014/11/24 11:31:30 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Generate/IsTwin.f90,v $  
!======================================================================!

  check=0

  do n=1,TwinN(n1,0)
    if(TwinN(n1,n) == n2) check=check+1
  end do

  do n=1,TwinN(n2,0)
    if(TwinN(n2,n) == n1) check=check+1
  end do

  if(check == 2) then
    IsTwin = .TRUE.
    return
  else if(check == 0) then
    IsTwin = .FALSE.
    return
  else
    write(*,*) 'IsTwin:   Major trouble !    Stopping !'
    stop
  endif

  END FUNCTION IsTwin
