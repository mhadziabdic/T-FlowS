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

!==============================================================================!
  subroutine Allocate_Unknown(phi, n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
!   This is to allocate an uknown for a normal time-integration algorithm.     !
!   Unknowns such as velocities and pressures should be allocated with it.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Unknown) :: phi
  integer       :: n_bnd_cells
  integer       :: n_cells
!------------------------------------------------------------------------------!

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-n_bnd_cells: n_cells));   phi % n  = 0.
  allocate (phi % o (-n_bnd_cells: n_cells));   phi % o  = 0.
  allocate (phi % oo(-n_bnd_cells: n_cells));   phi % oo = 0.

  ! Advection terms
  allocate (phi % C  (n_cells));   phi % C   = 0.
  allocate (phi % Co (n_cells));   phi % Co  = 0.
  allocate (phi % Coo(n_cells));   phi % coo = 0.

  ! Diffusion terms
  allocate (phi % Do (n_cells));   phi % Do  = 0.
  allocate (phi % Doo(n_cells));   phi % Doo = 0.

  ! Cross diffusion terms
  allocate (phi % X  (n_cells));   phi % X   = 0.
  allocate (phi % Xo (n_cells));   phi % Xo  = 0.
  allocate (phi % Xoo(n_cells));   phi % Xoo = 0.

  end subroutine

!==============================================================================!
  subroutine Allocate_Unknown_Statistics(phi, n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
!   This is to allocate additional values for statistics.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Unknown) :: phi
  integer       :: n_bnd_cells
  integer       :: n_cells
!------------------------------------------------------------------------------!

  ! Terms for statistics
  allocate (phi % mean(-n_bnd_cells: n_cells));   phi % mean = 0.
  allocate (phi % fluc(-n_bnd_cells: n_cells));   phi % fluc = 0.
  allocate (phi % filt(-n_bnd_cells: n_cells));   phi % filt = 0.

  end subroutine

!==============================================================================!
  subroutine Allocate_Unknown_Value_Simplified(phi, n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Unknown) :: phi
  integer       :: n_bnd_cells
  integer       :: n_cells
!------------------------------------------------------------------------------!

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-n_bnd_cells: n_cells));   phi % n  = 0.

  end subroutine

  end module 
