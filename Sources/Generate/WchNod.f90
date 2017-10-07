!======================================================================!
  integer function WchNod(c, n) 
!----------------------------------------------------------------------!
!   Returns the local number (1-8) of node n in cell c.                !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: c, n
!-------------------------------[Locals]-------------------------------!
  integer :: i, j
!======================================================================!

  WchNod=0
  if (c  < 0) then 
    write(*,*) 'Which node: Cell non existent !'
    return
  endif

!----- try the node only 
  do i=1,8
    if( CellN(c,i)  ==  n) then
      goto 10
    end if 
  end do     

!----- if it failed try his twins also
  do j=1,TwinN(n,0)
    do i=1,8
      if( CellN(c,i)  ==  TwinN(n,j)) then
        goto 10
      end if 
    end do     
  end do

  WchNod=0
  write(*,*) 'Which node: Trouble, node not found !'
  write(*,*) 'x, y, z = ', x_node(n), y_node(n), z_node(n)
  write(*,*) 'cell    = ', c, level(c)
  return

10   WchNod=i
  return

  end function WchNod
