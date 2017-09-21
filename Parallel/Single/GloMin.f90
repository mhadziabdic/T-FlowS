!======================================================================!
  SUBROUTINE GloMin(PHI) 
!----------------------------------------------------------------------!
!   Estimates global minimum among all processors.                     !
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
!  $Id: GloMin.f90,v 1.1 2014/11/24 11:39:07 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Parallel/Single/GloMin.f90,v $  
!======================================================================!

!================================================
      call MPI_ALLREDUCE      &               
!-----------------------------------+------------
	     (PHI,            & ! send buffer
	      PHInew,         & ! recv buffer 
	      1,              & ! length     
	      MPI_REAL,       & ! datatype  
	      MPI_MIN,        & ! operation 
	      MPI_COMM_WORLD, &             
	      error) 
!================================================

  PHI = PHInew

  END SUBROUTINE GloMin
