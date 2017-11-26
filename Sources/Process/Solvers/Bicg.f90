!==============================================================================!
  subroutine Bicg(A, x, r1,        &
                  prec,niter,tol,  &
                  ini_res,fin_res)
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. BiCG Method.          !
!------------------------------------------------------------------------------!
!   Allows preconditioning of the system by:                                   !
!     1. Diagonal preconditioning                                              !
!     2. Incomplete Cholesky preconditioning                                   !
!                                                                              !
!   The type of precondtioning is chosen by setting the variable prec to 0     !
!   (for no preconditioning), 1 (for diagonal preconditioning) or 2 (for       !
!   incomplete Cholesky preconditioning)                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use par_mod
  use Matrix_Mod
  use Solvers_Mod
  use Work_Mod, only: p1 => r_cell_01,  &
                      p2 => r_cell_02,  &
                      q1 => r_cell_03,  &
                      q2 => r_cell_04,  &
                      r2 => r_cell_05   
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Matrix_Type) :: A           
  real              :: x(-A % pnt_grid % n_bnd_cells : A % pnt_grid % n_cells)
  real              :: r1(A % pnt_grid % n_cells)    !  [A]{x}={r1}
  integer           :: prec,  niter       !  preconditioning
  real              :: tol                !  tolerance
  real              :: ini_res, fin_res   !  residual
!-----------------------------------[Locals]-----------------------------------!
  integer :: N, NB
  real    :: alfa, beta, rho, rhoold, bnrm2, error
  integer :: i, j, k, iter, sub
!==============================================================================!

  N  = A % pnt_grid % n_cells
  NB = A % pnt_grid % n_bnd_cells

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(N, A, D, prec) 

  !???????????????????????????????????!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !???????????????????????????????????!
  bnrm2=0.0
  do i=1,N
    bnrm2=bnrm2+r1(i)*r1(i)
  end do
  call glosum(bnrm2)  
  bnrm2=sqrt(bnrm2)

  if(bnrm2 < tol) then 
    iter=0
    goto 1
  end if  

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual(N, NB, A, x, r1) 

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  error=0.0
  do i=1,N
    error=error + r1(i)*r1(i)
  end do
  call glosum(error)
  error  = sqrt(error)  

  !---------------------------------------------------------------!
  !   Residual after the correction and before the new solution   !
  !---------------------------------------------------------------!
  ini_res=error 

  if(error < tol) then
    iter=0
    goto 1
  end if  

  !----------------------!
  !   Choose initial r   !
  !----------------------!
  do i=1,N
    r2(i)=r1(i)
  end do

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter=1,niter   

    !----------------------!  
    !    solve Mz  = r     !
    !    solve Mz = r      !
    !   (q instead of z)   !
    !----------------------!
    call Prec_Solve(N, NB, A, D, q1, r1(1), prec) 
    call Prec_Solve(N, NB, A, D, q2, r2(1), prec) 

    !-----------------!
    !   rho = (z,r)   !
    !-----------------!
    rho=0
    do i=1,N
      rho=rho+q1(i)*r2(i)
    end do
    call glosum(rho)

    if(iter == 1) then
      do i=1,N
        p1(i) = q1(i)
        p2(i) = q2(i)
      end do        
    else
      beta=rho/rhoold
      do i=1,N
        p1(i) = q1(i) + beta*p1(i)
        p2(i) = q2(i) + beta*p2(i)
      end do
    end if

    !-------------!
    !   q = A p   !
    !   q= A p    ! 
    !-------------!
    do i=1,N
      q1(i)  = 0.0                     
      q2(i) = 0.0                     
      do j=A % row(i), A % row(i+1)-1     
        k=A % col(j)                    
        q1(i) = q1(i) + A % val(j) * p1(k)   
        q2(i) = q2(i) + A % val(j) * p2(k)  
      end do
    end do
    call Exchange(A % pnt_grid, p1)
    call Exchange(A % pnt_grid, p2)
    do sub=1,n_proc
      if(NBBe(sub)  <=  NBBs(sub)) then
        do k=NBBs(sub),NBBe(sub),-1
          i=BufInd(k)
          q1(i) = q1(i) + A % bou(k)*p1(k)
          q2(i) = q2(i) + A % bou(k)*p2(k)
        end do
      end if
    end do

    !------------------------!
    !   alfa = (z,r)/(p,q)   !
    !------------------------!
    alfa=0.0
    do i=1,N
      alfa=alfa+p2(i)*q1(i)
    end do
    call glosum(alfa)
    alfa=rho/alfa

    !--------------------!
    !   x = x + alfa p   !
    !   r = r - alfa q   !
    !--------------------!
    do i=1,N
      x(i)  = x(i)  + alfa*p1(i)
      r1(i) = r1(i) - alfa*q1(i)
      r2(i) = r2(i) - alfa*q2(i)
    end do

    !???????????????????????!
    !   Check convergence   !
    !???????????????????????!
    error=0.0
    do i=1,N
      error=error+r1(i)*r1(i)
    end do  
    call glosum(error)
    error=sqrt(error)

    if(error < tol) goto 1

    rhoold=rho

  end do     ! iter

1 fin_res = error
  niter = iter

  end subroutine
