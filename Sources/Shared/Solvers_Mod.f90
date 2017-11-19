!==============================================================================!
  module Solvers_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Preconditioning "matrix" (D)
  type(Matrix_Type) :: D  ! preconditioning "matrix"

  ! Helping vectors and arrays for conjugate gradient type of solvers
  real, allocatable :: p1(:), p2(:)
  real, allocatable :: q1(:), q2(:), r2(:)
  real, allocatable :: u1(:), u2(:), v1(:), v2(:)
  real, allocatable :: u1_plus_q1(:)

  contains

  include 'Solvers_Mod/Allocate_Vectors.f90'

  end module 
