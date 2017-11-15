!==============================================================================!
  subroutine CalcShear(grid, ui, vi, wi, She)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the shear stress.                                !
!------------------------------------------------------------------------------!
!   She = sqrt( 2 * Sij * Sij )                                                !
!   Sij = 1/2 ( dUi/dXj + dUj/dXi )                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: ui(-grid % n_bnd_cells:grid % n_cells),  &
                     vi(-grid % n_bnd_cells:grid % n_cells),  &
                     wi(-grid % n_bnd_cells:grid % n_cells)
  real            :: She(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  call Exchange(grid, ui)
  call Exchange(grid, vi)
  call Exchange(grid, wi)

  !---------------!
  !   SGS terms   !
  !---------------!
  call GraPhi(grid, ui, 1, Ux, .TRUE.)  ! dU/dx
  call GraPhi(grid, ui, 2, Uy, .TRUE.)  ! dU/dy
  call GraPhi(grid, ui, 3, Uz, .TRUE.)  ! dU/dz
  call GraPhi(grid, vi, 1, Vx, .TRUE.)  ! dV/dx
  call GraPhi(grid, vi, 2, Vy, .TRUE.)  ! dV/dy
  call GraPhi(grid, vi, 3, Vz, .TRUE.)  ! dV/dz
  call GraPhi(grid, wi, 1, Wx, .TRUE.)  ! dW/dx
  call GraPhi(grid, wi, 2, Wy, .TRUE.)  ! dW/dy
  call GraPhi(grid, wi, 3, Wz, .TRUE.)  ! dW/dz

  do c = 1, grid % n_cells
    She(c) = Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) + &
             0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) + & 
             0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) + & 
             0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c)) 

     Vort(c) = - (0.5*(Vz(c) - Wy(c))*(Vz(c) - Wy(c)) + &
                  0.5*(Uz(c) - Wx(c))*(Uz(c) - Wx(c)) + &
                  0.5*(Vx(c) - Uy(c))*(Vx(c) - Uy(c)))

  end do 

  She  = sqrt(2.0 * She)
  Vort = sqrt(2.0 * abs(Vort))

  end subroutine
