!==============================================================================!
  subroutine User_Mod_Source(grid, phi, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
  type(Var_Type)     :: phi
  type(Matrix_Type)  :: a_matrix
  real, dimension(:) :: b_vector
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type)       :: might_be_useful_too
  real,    allocatable :: working_r(:), arrays_r(:), come_r(:), here_r(:)
  integer, allocatable :: working_i(:), arrays_i(:), come_i(:), here_i(:)
  real                 :: depends_r, on_r, the_r, case_r
  integer              :: depends_i, on_i, the_i, case_i
!------------------------------------------------------------------------------!

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name == 'T' ) then  

  end if
  
  !---------------------------------------------!
  !  Set source for turbulent kintetic energy   !
  !---------------------------------------------!
  if( phi % name == 'KIN' ) then  

  end if

  !---------------------------!
  !  Set source for epsilon   !
  !---------------------------!
  if( phi % name == 'EPS' ) then  

  end if

  end subroutine  ! fourth level comments
