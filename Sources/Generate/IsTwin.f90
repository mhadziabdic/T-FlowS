!======================================================================!
  logical function IsTwin(n1, n2) 
!----------------------------------------------------------------------!
!   Checks if the nodes are twins, i.e. are they shared on periodicity !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: n1, n2
!-------------------------------[Locals]-------------------------------!
  integer :: n, check
!======================================================================!

  check=0

  do n=1,TwinN(n1,0)
    if(TwinN(n1,n) == n2) check=check+1
  end do

  do n=1,TwinN(n2,0)
    if(TwinN(n2,n) == n1) check=check+1
  end do

  if(check == 2) then
    IsTwin = .TRUE.
    return
  else if(check == 0) then
    IsTwin = .FALSE.
    return
  else
    write(*,*) 'IsTwin:   Major trouble !    Stopping !'
    stop
  endif

  end function IsTwin
