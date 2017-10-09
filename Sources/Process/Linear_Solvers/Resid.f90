!======================================================================!
  subroutine Resid(N, NB, A, x, r1) 
!----------------------------------------------------------------------!
!   Calculates residuals.                                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use allt_mod, only: Matrix
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer      :: N, NB
  type(Matrix) :: A
  real         :: x(-NB:N), r1(N)             !  [A]{x}={r1}
!-------------------------------[Locals]-------------------------------!
  integer  :: i,j,k,sub
!======================================================================!

!+++++++++++++++++++++++++!
!     r = b - Ax          !
!     => Parallelized     ! 
!+++++++++++++++++++++++++!
  do i=1,N
    do j=A % col(i),A % col(i+1)-1     
      k = A % row(j)                 
      r1(i) = r1(i) - A % val(j) * x(k)  
    end do
  end do
!      call exchange(x) 
  do sub=1,n_proc
    if(NBBe(sub)  <=  NBBs(sub)) then
      do k=NBBs(sub),NBBe(sub),-1
        i=BufInd(k)
        r1(i) = r1(i) - A % bou(k)*x(k)
      end do
    end if
  end do

  end subroutine Resid
