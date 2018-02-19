!==============================================================================!
  subroutine Control_Mod_Heat_Transfer(verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: YES, NO
  use Flow_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('HEAT_TRANSFER', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val == 'YES' ) then
    heat_transfer = YES

  else if( val == 'NO' ) then
    heat_transfer = NO

  else
    print *, '# Unknown state for heat transfer: ', trim(val)
    print *, '# Exiting!'
    stop 

  end if

  end subroutine
