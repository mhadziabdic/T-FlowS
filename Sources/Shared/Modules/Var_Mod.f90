!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
  implicit none

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type         
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
    real              :: bound(128)          ! boundary values
    real              :: init(128)           ! initial values
    real              :: pro(11024)          ! inlfow profile
    real              :: Sigma               ! sigma
  end type Var_Type

  contains 

  include 'Var_Mod_Allocate_New.f90'
  include 'Var_Mod_Allocate_Solution.f90'
  include 'Var_Mod_Allocate_Statistics.f90'

  end module 
