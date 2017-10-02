!======================================================================!
  subroutine Prec1(N,NONZ,A,Acol,Arow,Adia,D,prec) 
!----------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CG Method.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer  :: N, NONZ      

  real     :: A(NONZ)
  integer  :: Acol(N),Adia(N)
  integer  :: Arow(NONZ)
  real     :: D(N) 

  integer  :: prec
!-------------------------------[Locals]-------------------------------!
  real     :: sum1
  integer  :: i, j, k
!======================================================================!
                 
!->>>
!      integer c 
!      do c=1,N
!        write(*,*) 'Cell: ', c
!        write(*,*) 'Width: ', Acol(c+1)-Acol(c)
!        write(*,'(3I7)') Acol(c), Adia(c), Acol(c+1)-1
!        write(*,*) 'Diag: ', A(Adia(c))
!        write(*,'(25F15.9)') ( A(j),     j=Acol(c),Acol(c+1)-1 )
!        write(*,'(25I7)') ( Arow(j), j=Acol(c),Acol(c+1)-1 )
!        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - -'
!      end do

!+++++++++++++++++++++++++!
!     preconditioning     !
!+++++++++++++++++++++++++!

!----- 1) diagonal preconditioning -----!
!     => Parallelized                   ! 
  if(prec == 1) then        
    do i=1,N                     
      D(i)=A(Adia(i))           
    end do                      

!----- 2) incomplete cholesky preconditioning -----!
  else if(prec == 2) then   
    do i=1,N
      sum1=A(Adia(i))          
      do j=Acol(i), Adia(i)-1         ! only lower traingular
        k=Arow(j)                    
        sum1= sum1- D(k) * A(j)*A(j)  
      end do
      D(i) = 1.0 / sum1                 ! BUG ?
    end do

!----- .) no preconditioning -----!
!     => Parallelized             ! 
  else                          
    do i=1,N
      D(i)=1.0
    end do
  end if 

  end subroutine Prec1
