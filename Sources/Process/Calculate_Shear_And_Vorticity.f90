!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(grid)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the shear stress.                                !
!------------------------------------------------------------------------------!
!   She = sqrt( 2 * Sij * Sij )                                                !
!   Sij = 1/2 ( dUi/dXj + dUj/dXi )                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use Comm_Mod
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
  
  call Comm_Mod_Exchange(grid, u % n)
  call Comm_Mod_Exchange(grid, v % n)
  call Comm_Mod_Exchange(grid, w % n)

  !---------------!
  !   SGS terms   !
  !---------------!
  call GraPhi(grid, u % n, 1, u % x, .true.)  ! dU/dx
  call GraPhi(grid, u % n, 2, u % y, .true.)  ! dU/dy
  call GraPhi(grid, u % n, 3, u % z, .true.)  ! dU/dz
  call GraPhi(grid, v % n, 1, v % x, .true.)  ! dV/dx
  call GraPhi(grid, v % n, 2, v % y, .true.)  ! dV/dy
  call GraPhi(grid, v % n, 3, v % z, .true.)  ! dV/dz
  call GraPhi(grid, w % n, 1, w % x, .true.)  ! dW/dx
  call GraPhi(grid, w % n, 2, w % y, .true.)  ! dW/dy
  call GraPhi(grid, w % n, 3, w % z, .true.)  ! dW/dz

  do c = 1, grid % n_cells
    shear(c) = u % x(c)**2                     &
             + v % y(c)**2                     &
             + w % z(c)**2                     &
             + 0.5 * (v % z(c) + w % y(c))**2  & 
             + 0.5 * (u % z(c) + w % x(c))**2  & 
             + 0.5 * (v % x(c) + u % y(c))**2

    vort(c) = - (  0.5*(v % z(c) - w % y(c))**2   &
                 + 0.5*(u % z(c) - w % x(c))**2   &
                 + 0.5*(v % x(c) - u % y(c))**2)

  end do 

  shear = sqrt(2.0 * shear)
  vort  = sqrt(2.0 * abs(vort))

  end subroutine
