!======================================================================!
  subroutine EndPar
!----------------------------------------------------------------------!
!   Ends parallel execution.                                           !
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Include]-------------------------------!
  include 'mpif.h'
!-------------------------------[Locals]-------------------------------!
  integer :: ERROR
!======================================================================!

  call MPI_FINALIZE(ERROR)

  end subroutine
