!==============================================================================!
  subroutine Info_Mod_Time_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use par_mod    
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer               :: i
  character(len=L_LINE) :: tmp
!==============================================================================!

  if (this_proc < 2) then

    print '(a83)', ' '
    print '(a83)', time_info % line_lead  

    ! Print only lines which have colon in the first column :-)
    do i=1,6
      print '(a83)', time_info % lines(i)
    end do

    print '(a83)', time_info % line_trail  
    print '(a83)', ' '

  end if
                 
  end subroutine
