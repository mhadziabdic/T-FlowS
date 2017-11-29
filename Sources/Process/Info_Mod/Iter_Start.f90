!==============================================================================!
  subroutine Info_Mod_Iter_Start()
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
! use pro_mod    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n  ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  ! Create frame
  do i = 1, L_LINE, L_LINE-1
    iter_info % line_lead (i:i) = '#'
    iter_info % lines(1)  (i:i) = '#'
    iter_info % lines(2)  (i:i) = '#'
    iter_info % lines(3)  (i:i) = '#'
    iter_info % lines(4)  (i:i) = '#'
    iter_info % line_trail(i:i) = '#'
  end do

  do i = 2, L_LINE-1
    iter_info % line_lead (i:i) = '='
    iter_info % line_trail(i:i) = '-'
  end do

  ! Create separators
  do i = 2, L_LINE-L_BOX, L_BOX
    write(iter_info % lines(1) (i:i+L_BOX-1), '(a20)') '|'
    write(iter_info % lines(2) (i:i+L_BOX-1), '(a20)') '|'
    write(iter_info % lines(3) (i:i+L_BOX-1), '(a20)') '|'
    write(iter_info % lines(4) (i:i+L_BOX-1), '(a20)') '|'
    write(iter_info % line_lead (i+L_BOX-1 :  &
                                 i+L_BOX-1),   '(a1)') '+'
    write(iter_info % line_trail(i+L_BOX-1 :  &
                                 i+L_BOX-1),   '(a1)') '+'
  end do

  end subroutine
