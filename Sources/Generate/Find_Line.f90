!==============================================================================!
  subroutine Find_Line(n1, n2, res) 
!------------------------------------------------------------------------------!
!   Searches for a smallest block where the line defined by n1-n2 is.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Domain_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)  :: n1, n2
  integer, intent(out) :: res
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, l1, l2
!==============================================================================!

  do b = 1, size(dom % blocks)
    do l1 = 1, 8
      do l2 = 1, 8
        if( (dom % blocks(b) % points(l1) == n1) .and. &
            (dom % blocks(b) % points(l2) == n2) ) then
          if( iabs(l2-l1) == 1 ) res = dom % blocks(b) % resolutions(1) 
          if( iabs(l2-l1) == 2 ) res = dom % blocks(b) % resolutions(2) 
          if( iabs(l2-l1) == 4 ) res = dom % blocks(b) % resolutions(3) 
          goto 1
        end if 
      end do
    end do     
  end do 

  write(*,*) '# Error message form Generator'
  write(*,*) '# You tried to define the line', n1, n2, ' but it'
  write(*,*) '# doesn''t exists in the block specifications.'
  write(*,*) '# Exiting !'
  stop

1 return

  end subroutine Find_Line
