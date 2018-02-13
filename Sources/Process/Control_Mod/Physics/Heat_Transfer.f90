!==============================================================================!
  subroutine Control_Mod_Heat_Transfer(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('HEAT_TRANSFER', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val.ne.'YES' .and. val.ne.'NO') then
    print *, '# Unknown state for heat transfer: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
