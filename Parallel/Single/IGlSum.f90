!======================================================================!
  SUBROUTINE IGlSum(PHI) 
!----------------------------------------------------------------------!
!   Estimates global summ among all processors.                        !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: PHI
!-------------------------------[Locals]-------------------------------!
  INTEGER :: PHInew
  INTEGER :: error
!--------------------------------[CVS]---------------------------------!
!  $Id: IGlSum.f90,v 1.1 2014/11/24 11:39:07 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Parallel/Single/IGlSum.f90,v $  
!======================================================================!

!================================================
      call MPI_ALLREDUCE      &               
!-----------------------------------+------------
	     (PHI,            & ! send buffer
	      PHInew,         & ! recv buffer 
	      1,              & ! length     
	      MPI_INTEGER,    & ! datatype  
	      MPI_SUM,        & ! operation 
	      MPI_COMM_WORLD, &             
	      error) 
!================================================

  PHI = PHInew

  END SUBROUTINE IGlSum
