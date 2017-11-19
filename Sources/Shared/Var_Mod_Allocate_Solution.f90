!==============================================================================!
  subroutine Var_Mod_Allocate_Solution(name_phi, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: name_phi
  type(Var_Type)   :: phi
  type(Grid_Type)  :: grid
!==============================================================================!
  
  ! Store variable name
  phi % name = name_phi

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-grid % n_bnd_cells: grid % n_cells));   phi % n  = 0.
  allocate (phi % o (-grid % n_bnd_cells: grid % n_cells));   phi % o  = 0.
  allocate (phi % oo(-grid % n_bnd_cells: grid % n_cells));   phi % oo = 0.

  ! Advection terms
  allocate (phi % C  (grid % n_cells));   phi % C   = 0.
  allocate (phi % Co (grid % n_cells));   phi % Co  = 0.
  allocate (phi % Coo(grid % n_cells));   phi % coo = 0.

  ! Diffusion terms
  allocate (phi % Do (grid % n_cells));   phi % Do  = 0.
  allocate (phi % Doo(grid % n_cells));   phi % Doo = 0.

  ! Cross diffusion terms
  allocate (phi % X  (grid % n_cells));   phi % X   = 0.
  allocate (phi % Xo (grid % n_cells));   phi % Xo  = 0.
  allocate (phi % Xoo(grid % n_cells));   phi % Xoo = 0.

  ! Variable's boundary flux
  allocate (phi % q(-grid % n_bnd_cells: -1));   phi % q  = 0.

  end subroutine
