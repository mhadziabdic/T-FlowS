!==============================================================================!
  subroutine Control_Mod_Time_Integration_Of_Innertia(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_OF_INNERTIA', 'linear',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val.ne.'LINEAR' .and. val.ne.'PARABOLIC') then
    print *, '# Unknown time-integration scheme for innertia: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
