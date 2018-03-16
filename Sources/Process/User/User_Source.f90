!==============================================================================!
  subroutine User_Source(grid, phi, a_matrix, b_vector)
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
  use Parameters_Mod
  use pro_mod
  use all_mod
  use Tokenizer_Mod
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
  integer              :: m, i, c
!------------------------------------------------------------------------------!

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 

  !---------------------------------------------------------------------!
  !  Correction of temperature equation with periodic boundary condition! 
  !---------------------------------------------------------------------!
  if( phi % name == 'T' ) then
    if(HOT==YES) then
      do m=1,grid % n_materials
        if(bulk(m) % p_drop_x > 1.0e-8.or.&
           bulk(m) % p_drop_y > 1.0e-8.or.&
           bulk(m) % p_drop_z > 1.0e-8.or.&
           bulk(m) % flux_x_o > 1.0e-8.or.&
           bulk(m) % flux_y_o > 1.0e-8.or.&
           bulk(m) % flux_z_o > 1.0e-8) then
          do c = 1, grid % n_cells 
            b(c)=b(c) - 1.0*Qflux * U % n(c)/bulk(1) % flux_x_o * grid % vol(c)
!            b(c)=b(c) -   2.0*3.14*Qflux * U % n(c)/bulk(1) % flux_x_o * grid % vol(c)
          end do
        end if
      end do
    end if
  end if 

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
