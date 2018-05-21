!==============================================================================!
  subroutine User_Mod_Allocate(grid)
!------------------------------------------------------------------------------!
!   This function allocates memory for user-defined variables.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer          :: n
  character(len=4) :: v_name
!==============================================================================!

  !--------------------------------------!
  !   Allocate memory for user scalars   !
  !--------------------------------------!
  allocate(user_scalar(n_user_scalars))

  !-------------------------------------!
  !   Browse through all user scalars   !
  !-------------------------------------!
  do n = 1, n_user_scalars

    ! Set variable name like this for the time being, later in control file
    v_name = 'C_00'
    write(v_name(3:4),'(i2.2)') n

    ! Allocate memory for user scalar
    call Var_Mod_Allocate_Solution(v_name, user_scalar(n), grid)

  end do

  !-------------------------------------!
  !   Allocate memory for user arrays   !
  !-------------------------------------!
  allocate(user_array(n_user_arrays, -grid % n_bnd_cells:grid % n_cells))
  do n = 1, n_user_arrays 
    user_array(n,:) = 0.
  end do

  end subroutine
