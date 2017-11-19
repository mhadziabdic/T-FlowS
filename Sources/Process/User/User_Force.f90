!==============================================================================!
  subroutine User_Force(grid, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for velocity.      !
!   It is called from "Compute_Velocity" function, just before calling the     !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
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
!------------------------------------------------------------------------------!

  !----------------------------------------------------! 
  !                                                    !
  !   Set source depending on the velocity component   !
  !                                                    !
  !----------------------------------------------------! 

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name == 'U' ) then  

  end if
  
  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name == 'V' ) then  

  end if
  
  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name == 'W' ) then  

  end if

  end subroutine  ! fourth level comments
