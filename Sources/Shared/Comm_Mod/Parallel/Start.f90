!==============================================================================!
  subroutine Comm_Mod_Start
!------------------------------------------------------------------------------!
!   Initializes parallel execution.                                            !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  CALL Mpi_Init(error)

  ! Get number of processors
  CALL Mpi_Comm_Size(MPI_COMM_WORLD, n_proc, error) 

  ! Get current processor number
  CALL Mpi_Comm_Rank(MPI_COMM_WORLD, this_proc, error)

  ! Use Fortran counting - from 1
  this_proc = this_proc + 1

  ! But if run is sequential, set the only processor to zero
  if(n_proc == 1) then
    n_proc = 0
    this_proc = 0
  endif

  end subroutine
