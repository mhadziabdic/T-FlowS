!==============================================================================!
  subroutine Prec_Solve(A, x, b, prec) 
!------------------------------------------------------------------------------!
! Solves the preconditioning system [D]{x}={b}                                 !
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
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Matrix_Type) :: A           
  real              :: x(-A % pnt_grid % n_bnd_cells : A % pnt_grid % n_cells)
  real              :: b( A % pnt_grid % n_cells)
  integer           :: prec
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, N, NB
  real    :: sum1
!==============================================================================!
           
  N  = A % pnt_grid % n_cells
  NB = A % pnt_grid % n_bnd_cells

  !---------------------------------! 
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec == 1) then
    do i=1,N
      x(i)=b(i)/D % val(D % dia(i))
    end do

  !--------------------------------------------! 
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec == 2) then

    ! Forward substitutionn
    do i=1,N
      sum1=b(i)
      do j=A % row(i),A % dia(i)-1  ! only the lower triangular
        k = A % col(j)             
        sum1 = sum1- A % val(j)*x(k)  
      end do
      x(i) = sum1 * D % val(D % dia(i))         ! BUG ?
    end do

    do i=1,N
      x(i) = x(i) / ( D % val(D % dia(i)) + TINY )
    end do

    ! Backward substitution
    do i=N,1,-1
      sum1=x(i)
      do j = A % dia(i)+1, A % row(i+1)-1 ! upper triangular 
        k = A % col(j)                  
        sum1 = sum1 - A % val(j)*x(k)      
      end do
      x(i) = sum1* D % val(D % dia(i))               ! BUG ?
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    do i=1,N
      x(i)=b(i)
    end do
  end if

  end subroutine
