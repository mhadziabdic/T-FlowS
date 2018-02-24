!==============================================================================!
  subroutine Calculate_Sgs_Hybrid(grid)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use allp_mod
  use Flow_Mod
  use Comm_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
! near(c) is the number of corresponding cell on the nearest wall.
! In case that, in parallel executions, the subdomain does not have 
! any nearwall cells, the near(c) is zero.
! near(c) is calculated in NearWallCells.f90, only ones in the beginig
! of a simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2 
  real    :: Nx, Ny, Nz
  real    :: Cs, R
  real    :: Stot, lf, UtauL, Uff 
  real    :: Utot, Unor, Utan, Apow, Bpow, nu, dely, yPlus 
  real    :: fun
!==============================================================================!

  do c = 1, grid % n_cells
    lf = grid % vol(c) ** ONE_THIRD    
    vis_t_sgs(c) = density                &
                 * (lf*lf)                &          ! delta^2 
                 * c_dyn(c)               &          ! c_dynamic   
                 * shear(c)      
  end do

  call Comm_Mod_Exchange(grid, vis_t_sgs)

  end subroutine
