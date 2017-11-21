!==============================================================================!
  subroutine Compute_Shear_And_Vorticity(grid)
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
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!
  
  call Exchange(grid, u % n)
  call Exchange(grid, v % n)
  call Exchange(grid, w % n)

  !---------------!
  !   SGS terms   !
  !---------------!
  call GraPhi(grid, u % n, 1, Ux, .TRUE.)  ! dU/dx
  call GraPhi(grid, u % n, 2, Uy, .TRUE.)  ! dU/dy
  call GraPhi(grid, u % n, 3, Uz, .TRUE.)  ! dU/dz
  call GraPhi(grid, v % n, 1, Vx, .TRUE.)  ! dV/dx
  call GraPhi(grid, v % n, 2, Vy, .TRUE.)  ! dV/dy
  call GraPhi(grid, v % n, 3, Vz, .TRUE.)  ! dV/dz
  call GraPhi(grid, w % n, 1, Wx, .TRUE.)  ! dW/dx
  call GraPhi(grid, w % n, 2, Wy, .TRUE.)  ! dW/dy
  call GraPhi(grid, w % n, 3, Wz, .TRUE.)  ! dW/dz

  do c = 1, grid % n_cells
    Shear(c) = Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) +  &
               0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) +      & 
               0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) +      & 
               0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c)) 

    Vort(c) = - (0.5*(Vz(c) - Wy(c))*(Vz(c) - Wy(c)) +   &
                 0.5*(Uz(c) - Wx(c))*(Uz(c) - Wx(c)) +   &
                 0.5*(Vx(c) - Uy(c))*(Vx(c) - Uy(c)))

  end do 

  Shear = sqrt(2.0 * Shear)
  Vort  = sqrt(2.0 * abs(Vort))

  end subroutine
