!==============================================================================!
!   Module for global variable definitions.  Variables defined here are        !
!   supposed to be accessible to all routines in all programs.                 !
!==============================================================================!
module all_mod

  use allp_mod

  implicit none

  !----------------------------------------!
  !   Variables for ease of input/output   !
  !----------------------------------------!
  character(len=80)  :: problem_name

  !-------------------------------------------!
  !   Logical quantities desribing the grid   !
  !-------------------------------------------!
  integer, allocatable :: material(:)     ! material markers

  integer, allocatable :: TypeBC(:)       ! type of boundary condition

end module 
