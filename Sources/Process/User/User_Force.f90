!==============================================================================!
  subroutine User_Force(grid, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for velocity.      !
!   It is called from "Compute_Velocity" function, just before calling the     !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod
  use Info_Mod
  use Parameters_Mod
  use Grid_Mod
  use Var_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
  type(Var_Type)     :: ui        ! velocity component
  type(Matrix_Type)  :: a_matrix  ! system matrix
  real, dimension(:) :: b_vector  ! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type)       :: might_be_useful_too
  real,    allocatable :: working_r(:), arrays_r(:), come_r(:), here_r(:)
  integer, allocatable :: working_i(:), arrays_i(:), come_i(:), here_i(:)
  real                 :: depends_r, on_r, the_r, case_r
  integer              :: depends_i, on_i, the_i, case_i
  integer              :: c
  real                 :: grav
!------------------------------------------------------------------------------!


  !----------------------------------------------------! 
  !                                                    !
  !   Set source depending on the velocity component   !
  !                                                    !
  !----------------------------------------------------! 

  if(BUOY == YES) then
  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
    if( ui % name == 'U' ) then  
      grav = grav_x      
    end if
  
  !-------------------------------------------------------!
  !  Set source for velocity component in "y" direction   !
  !-------------------------------------------------------!
    if( ui % name == 'V' ) then  
      grav = grav_y      
    end if
  
  !-------------------------------------------------------!
  !  Set source for velocity component in "z" direction   !
  !-------------------------------------------------------!
    if( ui % name == 'W' ) then  
      grav = grav_z      
    end if

    do c=1,grid % n_cells
      b(c) = b(c) - grav*(T%n(c) - Tref)*grid % vol(c) 
    end do
  end if

  end subroutine  ! fourth level comments
