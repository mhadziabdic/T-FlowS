subroutine Cgns_Mod_Get_Arrays_Dimensions_Par(idx, NN_or_NC)
!   Fetches correct dimensions for arrays in CGNS lib dependent functions      !
!------------------------------------------------------------------------------!
!   Arrays structure in CGNS parallel functions are strictly followings:       !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  !
!------------------------------------------------------------------------------!
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!   Connections:  |-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|...!
!----------------------------------[Modules]-----------------------------------!
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
  include 'mpif.h'
!-----------------------------------[Locals]-----------------------------------!
  integer :: idx, NN_or_NC
  integer :: Array_at_root(1:n_proc)
  integer :: tmp(1:n_proc)
  integer :: i, ier
!------------------------------------------------------------------------------!

  Array_at_root = 0

  call wait
  call mpi_gather(  &
    NN_or_NC,       & ! send number of nodes/cells, which is
    1,              & ! 1 element
    MPI_INTEGER8,   & ! of 64-bit integer type
    Array_at_root,  & ! to the root array
    1,              & ! which is 1 per processor long message.
    MPI_INTEGER8,   & ! Received data type is 64-bit integer
    0,              & ! at root processor
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error

  if (this_proc == 1) then

    tmp(:) = Array_at_root

    do i = 2, n_proc
      Array_at_root(i) = sum(tmp(1: i-1))
    end do
    Array_at_root(1) = 0
    Array_at_root = 1 + Array_at_root

  end if

  call wait
  call mpi_scatter( &
    Array_at_root,  & ! send number of nodes/cells, which is
    1,              & ! 1 elements per processor
    MPI_INTEGER8,   & ! of 64-bit integer type
    idx,            & ! to integer "idx"
    1,              & ! which is 1 element long message.
    MPI_INTEGER8,   & ! Received data type is 64-bit integer
    0,              & ! from the proc 1
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error

end subroutine