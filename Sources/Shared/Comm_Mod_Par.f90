!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
!==============================================================================!

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! These names are ugly but mean number of buffer boundaries start and end
  integer, allocatable :: nbb_s(:), nbb_e(:)

  integer, allocatable :: buffer_index(:)

  contains

  include 'Comm_Mod/Parallel/Allocate_Memory.f90'
  include 'Comm_Mod/Parallel/End.f90'
  include 'Comm_Mod/Parallel/Exchange.f90'
  include 'Comm_Mod/Parallel/Global_Max_Real.f90'
  include 'Comm_Mod/Parallel/Global_Min_Real.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Real.f90'
  include 'Comm_Mod/Parallel/Load_Buffers.f90'
  include 'Comm_Mod/Parallel/Start.f90'
  include 'Comm_Mod/Parallel/Wait.f90'

  end module
