!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!      for 'LES' computations       !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module les_mod

  use Var_Mod
  use Turbulence_Mod

  implicit none 

  ! Variables relevant for 'LES' computations
  real              :: ReTau, Cs0, Kflow  
  real, allocatable :: Utau(:), Vtau(:), Wtau(:)
  real, allocatable :: c_dyn(:), c_dyn_mean(:)

  ! Used in Dynamic Smgaorinsky model 
  real,allocatable :: Aval_dif(:)

  ! For LES you need to know nearest wall cell
  integer, allocatable :: nearest_wall_cell(:)

  real,allocatable :: shear_mean(:), kin_sgs(:)
  real,allocatable :: vis_t_sgs(:), vis_t_mean(:)

  real,allocatable :: shear_r(:), shear_mean_r(:), wale_v(:)
  real,allocatable ::                             &
                       UUf(:), VVf(:), WWf(:),    &
                       UVf(:), UWf(:), VWf(:),    &
                       M11f(:), M22f(:), M33f(:), &
                       M12f(:), M13f(:), M23f(:), &
                       shear_test(:),             &
                       Cinst(:)

end module  
