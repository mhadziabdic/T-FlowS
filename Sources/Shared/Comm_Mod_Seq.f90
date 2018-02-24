!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for no MPI functionality.                                           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! These names are ugly but mean number of buffer boundaries start and end
  integer, allocatable :: nbb_s(:), nbb_e(:)

  integer, allocatable :: BufInd(:)

  contains

  include 'Comm_Mod/Sequential/End.f90'
  include 'Comm_Mod/Sequential/Exchange.f90'
  include 'Comm_Mod/Sequential/Global_Max_Real.f90'
  include 'Comm_Mod/Sequential/Global_Min_Real.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real.f90'
  include 'Comm_Mod/Sequential/Load_Buffers.f90'
  include 'Comm_Mod/Sequential/Start.f90'
  include 'Comm_Mod/Sequential/Wait.f90'

  end module
