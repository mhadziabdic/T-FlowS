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

    print *, bulk_info % line_lead  
    print *, bulk_info % lines(1)
    print *, bulk_info % line_sep
    print *, bulk_info % lines(2)
    print *, bulk_info % lines(3)
    print *, bulk_info % line_trail  
    print *, ' '

  end if
                 
  end subroutine
