!==============================================================================!
  subroutine Find_Surface(n1, n2, n3, n4, block, face) 
!------------------------------------------------------------------------------!
!   Searches for a block where the surface defined by n1, n2, n3, n4 is.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)  :: n1, n2, n3, n4
  integer, intent(out) :: block, face
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, fc, p1, p2, p3, p4
!==============================================================================!

  do b=1,Nbloc
    do fc=1,6
      p1=block_faces(b, fc, 1)
      p2=block_faces(b, fc, 2)
      p3=block_faces(b, fc, 3)
      p4=block_faces(b, fc, 4) 
      if( ((p1 == n1).and.(p3 == n3)) .or.  &
          ((p1 == n4).and.(p3 == n2)) .or.  &
          ((p1 == n3).and.(p3 == n1)) .or.  &
          ((p1 == n2).and.(p3 == n4)) ) goto 1
    end do     
  end do 

  write(*,*) '# Error message from Generator'
  write(*,*) '# You tried to define the surface', n1, n2, n3, n4
  write(*,*) '# but it doesn''t exists in the block specifications.'
  write(*,*) '# Exiting !'
  stop

1 block=b
  face =fc

  end subroutine Find_Surface
