!==============================================================================!
  subroutine User_Mod_Save_Results(grid, current_time_step)
!------------------------------------------------------------------------------!
!   This is a prototype of customized saving functions.  All cases which       !
!   require special procedures when saving results (1D or 2D profiles, non-    !
!   standard save formats, and alike) should be done here, case-based.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: current_time_step
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type)       :: might_be_useful_too
  real,    allocatable :: working_r(:), arrays_r(:), come_r(:), here_r(:)
  integer, allocatable :: working_i(:), arrays_i(:), come_i(:), here_i(:)
  real                 :: depends_r, on_r, the_r, case_r
  integer              :: depends_i, on_i, the_i, case_i
!------------------------------------------------------------------------------!

  !--------------------------! 
  !   First level comments   !
  !--------------------------! 

  !---------------------------! 
  !   Second level comments   !
  !---------------------------! 

  ! Third level comments

  end subroutine  ! fourth level comments
