!==============================================================================!
  subroutine Info_Mod_Bulk_Print()
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

    print '(a83)', bulk_info % line_lead  
    print '(a83)', bulk_info % lines(1)
    print '(a83)', bulk_info % line_sep
    print '(a83)', bulk_info % lines(2)
    print '(a83)', bulk_info % lines(3)
    print '(a83)', bulk_info % line_trail  
    print '(a83)', ' '

  end if
                 
  end subroutine
