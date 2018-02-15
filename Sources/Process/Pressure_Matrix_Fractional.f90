!==============================================================================!
  subroutine Pressure_Matrix_Fractional(grid, dt)
!------------------------------------------------------------------------------!
!   Forms the pressure system matrix for the fractional step method.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
!-----------------------------------[Locals]-----------------------------------!
  real    :: A12
  integer :: c, c1, c2, s 
!==============================================================================!

  do c = 1, A % row(grid % n_cells+1)  ! this is number of nozero entries + 1
    A % val(c) = 0.0
  end do

  !-----------------------------!
  !   Calculate system matrix   ! 
  !-----------------------------!
  do s = 1, grid % n_faces    

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  > 0) then
      A12 = dt * Scoef(s) 
      A % val(A % pos(1,s)) = -A12
      A % val(A % pos(2,s)) = -A12
      A % val(A % dia(c1)) =                                            &
      A % val(A % dia(c1)) +  A12
      A % val(A % dia(c2)) =                                            &
      A % val(A % dia(c2)) +  A12
    else
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then
        A12 = dt * Scoef(s)
        A % val(A % dia(c1)) =                                          &
        A % val(A % dia(c1)) +  A12
        A % bou(c2) = -A12
      end if
    end if 

  end do ! through faces

  end subroutine
