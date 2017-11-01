!==============================================================================!
  logical function Are_Nodes_Twins(n1, n2) 
!------------------------------------------------------------------------------!
!   Checks if the nodes are twins, i.e. are they shared on periodicity         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod, only: TwinN
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n1, n2
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, check
!==============================================================================!

  check=0

  do n=1,TwinN(n1,0)
    if(TwinN(n1,n) == n2) check=check+1
  end do

  do n=1,TwinN(n2,0)
    if(TwinN(n2,n) == n1) check=check+1
  end do

  if(check == 2) then
    Are_Nodes_Twins = .TRUE.
    return
  else if(check == 0) then
    Are_Nodes_Twins = .FALSE.
    return
  else
    write(*,*) 'Are_Nodes_Twins:   Major trouble !    Stopping !'
    stop
  endif

  end function Are_Nodes_Twins
