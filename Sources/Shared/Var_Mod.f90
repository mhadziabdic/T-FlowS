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
    character(len=4)  :: name                ! variable name, always make it 
                                             ! uppercase and keep very short
    real, allocatable :: n(:)                ! new value
    real, allocatable :: o(:), oo(:)         ! old and older then old
    real, allocatable :: C(:), Co(:), Coo(:) ! convective fluxes
    real, allocatable :: Do(:), Doo(:)       ! difussive fluxes
    real, allocatable :: X(:), Xo(:), Xoo(:) ! surfce sources  
    real, allocatable :: mean(:)             ! long time average
    real, allocatable :: filt(:)             ! long time average
    real, allocatable :: fluc(:) 
    real, allocatable :: q(:)                ! flux of a variable
    real              :: URF                 ! under relaxation factor
    real              :: Stol                ! solver tolerance
    real              :: bound(1024)         ! boundary values
    real              :: init(1024)          ! initial values
    real              :: pro(11024)          ! inlfow profile
    real              :: Sigma               ! sigma
    type(Grid_Type), pointer :: pnt_grid
  end type

  contains 

  include 'Var_Mod_Allocate_New_Only.f90'
  include 'Var_Mod_Allocate_Solution.f90'
  include 'Var_Mod_Allocate_Statistics.f90'

  end module 
