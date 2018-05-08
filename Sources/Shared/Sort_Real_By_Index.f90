!======================================================================!
  subroutine Sort_Real_By_Index(X,indx,N)
!----------------------------------------------------------------------!
!   Sorts real array X according to indx.                              !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer :: N,indx(N)
  real    :: X(N)
!-------------------------------[Locals]-------------------------------!
  integer           :: i
  real, allocatable :: work(:)
!======================================================================!

  allocate(work(N)); work=0.

  do i=1,N
    work(indx(i))=X(i)
  end do

  do i=1,N
    X(i)=work(i)
  end do

  deallocate(work)

  end subroutine
