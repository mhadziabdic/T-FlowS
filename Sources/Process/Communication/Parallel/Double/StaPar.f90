!======================================================================!
  subroutine StaPar
!----------------------------------------------------------------------!
!   Initializes parallel execution.                                    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use par_mod 
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Include]-------------------------------!
  include 'mpif.h'
!-------------------------------[Locals]-------------------------------!
  integer :: error
!======================================================================!

  write(*,*) 'Calling MPI_INIT'

  CALL MPI_INIT(ERROR)

  write(*,*) 'Calling MPI_COMM_WORLD'

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,                                &
                     n_proc,                                          &
                     error) 

  write(*,*) 'Calling MPI_COMM_RANK'

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,                                &
                     this_proc,                                          &
                     error)

  this_proc=this_proc+1

  if(n_proc == 1) then
    n_proc = 0
    this_proc = 0
  endif

  end subroutine
