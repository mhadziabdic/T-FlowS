!==============================================================================!
  subroutine Info_Mod_Time_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer               :: i
  character(len=L_LINE) :: tmp
!==============================================================================!

  print *, ' '
  print *, time_info % line_lead  

  ! Print only lines which have colon in the first column :-)
  do i=1,6
    print *, time_info % lines(i)
  end do

  print *, time_info % line_trail  
  print *, ' '
                 
  end subroutine
