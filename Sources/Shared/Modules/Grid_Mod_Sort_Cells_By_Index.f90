!==============================================================================!
  subroutine Sort_Cells_By_Index(cells, indx, n)
!------------------------------------------------------------------------------!
!   Sorts array of cells according to indx.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Cell_Type) :: cells(n)
  integer         :: n, indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer                      :: i
  type(Cell_Type), allocatable :: work(:)
!==============================================================================!

  allocate(work(n))

  do i = 1, n
    work( indx(i) ) % n_nodes = cells(i) % n_nodes
    work( indx(i) ) % n       = cells(i) % n
    work( indx(i) ) % c       = cells(i) % c
  end do

  do i=1,N
    cells(i) % n_nodes = work(i) % n_nodes
    cells(i) % n       = work(i) % n      
    cells(i) % c       = work(i) % c        
  end do

  deallocate(work)

  end subroutine
