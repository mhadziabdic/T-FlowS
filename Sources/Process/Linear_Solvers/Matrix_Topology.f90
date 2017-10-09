!======================================================================!
  subroutine Matrix_Topology(M)
!----------------------------------------------------------------------!
!   Determines the topology of the system matrix.                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use allt_mod, only: Matrix
  use all_mod,  only: NC, NS, SideC
  use sol_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  type(Matrix) :: M
!-------------------------------[Locals]-------------------------------!
  integer              :: c, s, j, n
  integer              :: c1, c2
  integer, allocatable :: stencw(:)
!======================================================================!
!   Relies only on SideC structure. Try to keep it that way.
!----------------------------------------------------------------------!
                  
!---- memory allocation
  allocate(stencw(NC)); stencw=1

!====================================================================!
!                                                                    !
!     Determine the topology of unstructured system of equations     !
!                                                                    !
!====================================================================!

  if(this_proc < 2) write(*,*) '# Determining matrix topology.'

  ! Initialize number of nonzeros
  M % nonzeros = 0

  ! Compute stencis widths
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
      stencw(c1)=stencw(c1)+1
      stencw(c2)=stencw(c2)+1
    end if      
  end do

  ! Count the nonzero entries and allocate the memory for the array 
  n = 0  
  do c=1,NC
    n = n + stencw(c)
  end do   
  M % nonzeros = n + 1
  allocate(M % val(n+1)); M % val=0 ! it reffers to M % col+1 
  allocate(M % row(n+1)); M % row=0 ! it reffers to M % col+1 

  ! Form M % col and diagonal only formation of M % row
  M % col(1)=1
  do c=1,NC
    M % col(c+1)=M % col(c)+stencw(c)
    M % row(M % col(c)) = c   ! sam sebi je prvi
    stencw(c)=1
  end do 

  ! Extend M % row entries with off-diagonal terms      
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
      M % row(M % col(c1)+stencw(c1)) = c2
      M % row(M % col(c2)+stencw(c2)) = c1
      stencw(c1)=stencw(c1)+1
      stencw(c2)=stencw(c2)+1
    end if      
  end do

  ! Sort M % row to make them nice and neat
  do c=1,NC
    call isort(M % row(M % col(c)),                                       &
               M % row(M % col(c)),stencw(c),1)
    do j=M % col(c),M % col(c+1)-1
      if(M % row(j) == c) M % dia(c)=j
    end do
  end do 

  ! Connect faces with matrix entries
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then

      ! Where is M(c1,c2) and ...
      do c=M % col(c1),M % col(c1+1)-1 
        if(M % row(c)  ==  c2) M % pos(1,s)=c
      end do

      ! Where is M(c2,c1) and ...
      do c=M % col(c2),M % col(c2+1)-1 
        if(M % row(c)  ==  c1) M % pos(2,s)=c
      end do
    end if
  end do

  if(this_proc < 2) write(*,*) '# Finished !'
 
  deallocate(stencw)

  RETURN

  end subroutine Matrix_Topology
