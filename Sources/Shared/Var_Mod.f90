!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type         
    character(len=4)  :: name                  ! variable name, always make it 
                                               ! uppercase and keep very short
    real, allocatable :: n(:)                  ! new value
    real, allocatable :: o(:), oo(:)           ! old and older then old
    real, allocatable :: a(:), a_o(:), a_oo(:) ! advection fluxes
    real, allocatable :: d_o(:), d_oo(:)       ! difussion fluxes
    real, allocatable :: c(:), c_o(:), c_oo(:) ! cross-difusion
    real, allocatable :: mean(:)               ! time average
    real, allocatable :: filt(:)               ! filtered quantity
    real, allocatable :: fluc(:)               ! fluctuating value
    real, allocatable :: x(:), y(:), z(:)      ! gradient components
    real, allocatable :: q(:)                  ! flux of a variable
    real              :: URF                   ! under relaxation factor
    real              :: Stol                  ! solver tolerance
    real              :: bound(1024)           ! boundary values
    real              :: init(1024)            ! initial values
    real              :: pro(11024)            ! inlfow profile
    real              :: Sigma                 ! sigma
    type(Grid_Type), pointer :: pnt_grid
  end type

  contains 

  include 'Var_Mod/Allocate_New_Only.f90'
  include 'Var_Mod/Allocate_Solution.f90'
  include 'Var_Mod/Allocate_Statistics.f90'
  include 'Var_Mod/Allocate_Gradients.f90'

  end module 
