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
!-----------------------------[Arguments]------------------------------!
  type(Matrix) :: M
!-------------------------------[Locals]-------------------------------!
  integer              :: c, s, j, n
  integer              :: c1, c2
  integer, allocatable :: stencw(:)
!======================================================================!
!   Relies only on SideC structure. Try to keep it that way.
!----------------------------------------------------------------------!
                  
  ! Memory allocation
  allocate(stencw(NC)); stencw=1

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
  allocate(M % val(n+1)); M % val=0 ! it reffers to M % row+1 
  allocate(M % col(n+1)); M % col=0 ! it reffers to M % row+1 

  ! Form M % row and diagonal only formation of M % col
  M % row(1)=1
  do c=1,NC
    M % row(c+1)=M % row(c)+stencw(c)
    M % col(M % row(c)) = c   ! sam sebi je prvi
    stencw(c)=1
  end do 

  ! Extend M % col entries with off-diagonal terms      
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then
      M % col(M % row(c1)+stencw(c1)) = c2
      M % col(M % row(c2)+stencw(c2)) = c1
      stencw(c1)=stencw(c1)+1
      stencw(c2)=stencw(c2)+1
    end if      
  end do

  ! Sort M % col to make them nice and neat
  do c=1,NC
    call Sort_Int_Carry_Int(M % col(M % row(c)),  &
               M % col(M % row(c)),stencw(c),1)
    do j=M % row(c),M % row(c+1)-1
      if(M % col(j) == c) M % dia(c)=j
    end do
  end do 

  ! Connect faces with matrix entries
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  > 0) then

      ! Where is M(c1,c2) and ...
      do c=M % row(c1),M % row(c1+1)-1 
        if(M % col(c)  ==  c2) M % pos(1,s)=c
      end do

      ! Where is M(c2,c1) and ...
      do c=M % row(c2),M % row(c2+1)-1 
        if(M % col(c)  ==  c1) M % pos(2,s)=c
      end do
    end if
  end do

  if(this_proc < 2) write(*,*) '# Finished !'
 
  deallocate(stencw)

  RETURN

  end subroutine Matrix_Topology
