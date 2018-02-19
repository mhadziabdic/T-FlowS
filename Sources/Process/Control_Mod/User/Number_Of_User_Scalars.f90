!==============================================================================!
  subroutine Control_Mod_Number_Of_User_Scalars(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_USER_SCALARS', 0, &
                                  val, verbose)

  end subroutine
