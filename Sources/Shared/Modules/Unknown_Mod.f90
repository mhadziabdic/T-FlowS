!==============================================================================!
  module Unknown_Mod
!------------------------------------------------------------------------------!
  implicit none

  !------------------!
  !   Unknown type   !
  !------------------!
  type Unknown          
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
  end type Unknown

  contains 

  include 'Unknown_Mod_Allocate_Unknown.f90'
  include 'Unknown_Mod_Allocate_Unknown_New_Value_Only.f90'
  include 'Unknown_Mod_Allocate_Unknown_Statistics.f90'

  end module 
