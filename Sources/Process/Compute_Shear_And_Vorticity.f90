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

  call GraPhi(grid, u % n, 1, u % x, .TRUE.)  ! dU/dx
  call GraPhi(grid, u % n, 2, u % y, .TRUE.)  ! dU/dy
  call GraPhi(grid, u % n, 3, u % z, .TRUE.)  ! dU/dz
  call GraPhi(grid, v % n, 1, v % x, .TRUE.)  ! dV/dx
  call GraPhi(grid, v % n, 2, v % y, .TRUE.)  ! dV/dy
  call GraPhi(grid, v % n, 3, v % z, .TRUE.)  ! dV/dz
  call GraPhi(grid, w % n, 1, w % x, .TRUE.)  ! dW/dx
  call GraPhi(grid, w % n, 2, w % y, .TRUE.)  ! dW/dy
  call GraPhi(grid, w % n, 3, w % z, .TRUE.)  ! dW/dz

  do c = 1, grid % n_cells
    Shear(c) = u % x(c)**2                     &
             + v % y(c)**2                     &
             + w % z(c)**2                     &
             + 0.5 * (v % z(c) + w % y(c))**2  & 
             + 0.5 * (u % z(c) + w % x(c))**2  & 
             + 0.5 * (v % x(c) + u % y(c))**2

    Vort(c) = - (  0.5*(v % z(c) - w % y(c))**2   &
                 + 0.5*(u % z(c) - w % x(c))**2   &
                 + 0.5*(v % x(c) - u % y(c))**2)

  end do 

  Shear = sqrt(2.0 * Shear)
  Vort  = sqrt(2.0 * abs(Vort))

  end subroutine
