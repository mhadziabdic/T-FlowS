!======================================================================!
  subroutine INSort(X,indx,N)
!----------------------------------------------------------------------!
!   Sorts int. array X according to indx.                              !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer :: N,X(N),indx(N)
!-------------------------------[Locals]-------------------------------!
  integer             :: i
  integer,allocatable :: work(:)
!======================================================================!

  allocate(work(N)); work=0

  do i=1,N
    work(indx(i))=X(i)
  end do

  do i=1,N
    X(i)=work(i)
  end do

  deallocate(work)

  end subroutine INSort
