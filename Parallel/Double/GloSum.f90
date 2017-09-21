!======================================================================!
  SUBROUTINE GloSum(PHI) 
!----------------------------------------------------------------------!
!   Estimates global summ among all processors.                        !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  REAL    :: PHI
!-------------------------------[Locals]-------------------------------!
  REAL    :: PHInew
  INTEGER :: error
!--------------------------------[CVS]---------------------------------!
!  $Id: GloSum.f90,v 1.1 2014/11/24 11:39:27 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Parallel/Double/GloSum.f90,v $  
!======================================================================!

!================================================
      call MPI_ALLREDUCE      &               
!-----------------------------------+------------
	     (PHI,            & ! send buffer
	      PHInew,         & ! recv buffer 
	      1,              & ! length     
	      MPI_DOUBLE_PRECISION,     & ! datatype  
	      MPI_SUM,        & ! operation 
	      MPI_COMM_WORLD, &             
	      error) 
!================================================

  PHI = PHInew

  END SUBROUTINE GloSum
