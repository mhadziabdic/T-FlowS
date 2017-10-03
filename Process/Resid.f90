!======================================================================!
  subroutine Resid(N, NB, NONZ, A, Acol,Arow,Ab,x,r1) 
!----------------------------------------------------------------------!
!   Calculates residuals.                                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer  :: N, NB, NONZ      

  real     :: A(NONZ),Ab(-NB:-1)
  integer  :: Acol(N+1)
  integer  :: Arow(NONZ)
  real     :: x(-NB:N), r1(N)             !  [A]{x}={r1}
!-------------------------------[Locals]-------------------------------!
  integer  :: i,j,k,sub
!======================================================================!

!+++++++++++++++++++++++++!
!     r = b - Ax          !
!     => Parallelized     ! 
!+++++++++++++++++++++++++!
  do i=1,N
    do j=Acol(i),Acol(i+1)-1     
      k = Arow(j)                 
      r1(i) = r1(i) - A(j) * x(k)  
    end do
  end do
!      call exchange(x) 
  do sub=1,Npro
    if(NBBe(sub)  <=  NBBs(sub)) then
      do k=NBBs(sub),NBBe(sub),-1
        i=BufInd(k)
        r1(i) = r1(i) - Ab(k)*x(k)
      end do
    end if
  end do

  end subroutine Resid
