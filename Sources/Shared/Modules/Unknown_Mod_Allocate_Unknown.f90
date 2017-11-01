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
