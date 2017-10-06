MODULE sol_mod 

  use allp_mod

  implicit none

  real,allocatable :: D(:)
  real,allocatable :: p1(:), p2(:)
  real,allocatable :: q1(:), q2(:), r2(:)
  real,allocatable :: u1(:), u2(:), v1(:), v2(:)
  real :: alfa, beta, rho, rhoold, bnrm2, sum1, sum2, error

  integer :: i, j, k, iter, sub
end MODULE 
