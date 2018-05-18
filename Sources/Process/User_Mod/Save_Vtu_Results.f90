!==============================================================================!
  subroutine User_Mod_Save_Vtu_Results(grid)
!------------------------------------------------------------------------------!
!   This function is called from Save_Vtu and allows user to save his          !
!   own arrays in the Vtu file with results.                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
!----------------------------------[Locals]------------------------------------!
  integer          :: n
  character(len=4) :: a_name
!-----------------------------[Local parameters]-------------------------------!
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  !-----------------------!
  !   Save user scalars   !
  !-----------------------!
  do n = 1, n_user_scalars
    call Save_Vtu_Scalar(grid, IN_4, IN_5, user_scalar(n) % name,  &
                                           user_scalar(n) % n)
  end do
  
  !----------------------!
  !   Save user arrays   !
  !----------------------!
  do n = 1, n_user_arrays  

    a_name = 'UA00'
    write(a_name(3:4), '(i2.2)') n

    call Save_Vtu_Scalar(grid, IN_4, IN_5, a_name, user_array(n,:))
  end do


  end subroutine  ! fourth level comments
