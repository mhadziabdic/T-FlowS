!==============================================================================!
  subroutine Solvers_Mod_Allocate_Vectors(n_bnd_cells, n_cells)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n_bnd_cells, n_cells
!==============================================================================!

  ! Allocate memory for vectors used in Solvers_Mod
  allocate (p1(-n_bnd_cells:n_cells));       p1=0
  allocate (p2(-n_bnd_cells:n_cells));       p2=0
  allocate (q1(-n_bnd_cells:n_cells));       q1=0
  allocate (q2(-n_bnd_cells:n_cells));       q2=0
  allocate (r2(n_cells));                    r2=0
  allocate (u1(n_cells));                    u1=0
  allocate (u2(n_cells));                    u2=0
  allocate (v1(n_cells));                    v1=0
  allocate (v2(n_cells));                    v2=0
  allocate (u1_plus_q1(n_cells));            u1_plus_q1=0

  end subroutine
