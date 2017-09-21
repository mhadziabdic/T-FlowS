!======================================================================!
  SUBROUTINE StaPar
!----------------------------------------------------------------------!
!   Initializes parallel execution.                                    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE par_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-------------------------------[Locals]-------------------------------!
  INTEGER :: error
!--------------------------------[CVS]---------------------------------!
!  $Id: StaPar.f90,v 1.1 2014/11/24 11:39:27 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Parallel/Double/StaPar.f90,v $  
!======================================================================!

  CALL MPI_INIT(ERROR)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,                                &
		     NPro,                                          &
		     error) 

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,                                &
		     this,                                          &
		     error)

  this=this+1

  if(Npro == 1) then
    Npro = 0
    this = 0
  endif

  END SUBROUTINE StaPar
