!==============================================================================!
  subroutine CalcShear(grid, Ui, Vi, Wi, She)
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
  real            :: Ui(-NbC:NC), Vi(-NbC:NC), Wi(-NbC:NC)
  real            :: She(-NbC:NC)
!-----------------------------------[Locals]-------------------------------!
  integer :: c, i
  real    :: Sii, Sjk 
!==============================================================================!

  call Exchng(Ui)
  call Exchng(Vi)
  call Exchng(Wi)

  !---------------!
  !   SGS terms   !
  !---------------!
  call GraPhi(grid, Ui, 1, Ux, .TRUE.)  ! dU/dx
  call GraPhi(grid, Ui, 2, Uy, .TRUE.)  ! dU/dy
  call GraPhi(grid, Ui, 3, Uz, .TRUE.)  ! dU/dz
  call GraPhi(grid, Vi, 1, Vx, .TRUE.)  ! dV/dx
  call GraPhi(grid, Vi, 2, Vy, .TRUE.)  ! dV/dy
  call GraPhi(grid, Vi, 3, Vz, .TRUE.)  ! dV/dz
  call GraPhi(grid, Wi, 1, Wx, .TRUE.)  ! dW/dx
  call GraPhi(grid, Wi, 2, Wy, .TRUE.)  ! dW/dy
  call GraPhi(grid, Wi, 3, Wz, .TRUE.)  ! dW/dz

  do c=1,NC
    She(c) = Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) + &
             0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) + & 
             0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) + & 
             0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c)) 

     Vort(c) = - (0.5*(Vz(c) - Wy(c))*(Vz(c) - Wy(c)) + &
                  0.5*(Uz(c) - Wx(c))*(Uz(c) - Wx(c)) + &
                  0.5*(Vx(c) - Uy(c))*(Vx(c) - Uy(c)))

  end do 

  She = sqrt(2.0 * She)
  Vort = sqrt(2.0*abs(Vort))

  end subroutine
