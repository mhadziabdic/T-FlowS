!==============================================================================!
  subroutine Control_Mod_Save_Restart_Name(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SAVE_RESTART_NAME', 'skip',  &
                                   val, verbose)

  end subroutine
