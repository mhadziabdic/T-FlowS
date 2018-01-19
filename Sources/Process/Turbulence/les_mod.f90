!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!      for LES computations       !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module les_mod

  use Var_Mod
  use Turbulence_Mod

  implicit none 

  ! Variables relevant for LES computations
  real             :: ReTau, Cs0, Kflow  
  real,allocatable :: Utau(:), Vtau(:), Wtau(:)
  real,allocatable :: Cdyn(:), Cdyn_mean(:)

  ! Used in Dynamic Smgaorinsky model 
  real,allocatable :: Aval_dif(:)

  real,allocatable :: Shear(:), ShearMean(:), Ksgs(:), TauWall(:), VISt_mean(:)
  real,allocatable :: Shear_r(:), ShearMean_r(:), WALEv(:)
  real,allocatable ::                             &
                       UUf(:), VVf(:), WWf(:),    &
                       UVf(:), UWf(:), VWf(:),    &
                       M11f(:), M22f(:), M33f(:), &
                       M12f(:), M13f(:), M23f(:), &
                       ShearTest(:),              &
                       Cinst(:)
  real,allocatable ::  Puu_mean(:), Pvv_mean(:), Pww_mean(:), &
                       Puv_mean(:), Puw_mean(:), Pvw_mean(:), &
                       Put_mean(:), Pvt_mean(:), Pwt_mean(:), Ptt_mean(:), &
                       Diss_uu_mean(:), Diss_vv_mean(:), Diss_ww_mean(:), Diss_sgs_mean(:), &
                       Diss_uv_mean(:), Diss_uw_mean(:), Diss_vw_mean(:), &
                       Diss_ut_mean(:), Diss_vt_mean(:), Diss_wt_mean(:), Diss_tt_mean(:),&
                       Difv_uu_mean(:), Difv_vv_mean(:), Difv_ww_mean(:), &
                       Difv_uv_mean(:), Difv_uw_mean(:), Difv_vw_mean(:), &
                       Difv_ut_mean(:), Difv_vt_mean(:), Difv_wt_mean(:), Difv_tt_mean(:),&
                       Dift_uu_mean(:), Dift_vv_mean(:), Dift_ww_mean(:), &
                       Dift_uv_mean(:), Dift_uw_mean(:), Dift_vw_mean(:), &
                       Dift_ut_mean(:), Dift_vt_mean(:), Dift_wt_mean(:), Dift_tt_mean(:),&
                       PR_uu_mean(:), PR_vv_mean(:), PR_ww_mean(:), &
                       PR_uv_mean(:), PR_uw_mean(:), PR_vw_mean(:), &
                       PR_ut_mean(:), PR_vt_mean(:), PR_wt_mean(:), PR_tt_mean(:),&
                       PD_uu_mean(:), PD_vv_mean(:), PD_ww_mean(:), &
                       PD_uv_mean(:), PD_uw_mean(:), PD_vw_mean(:), &
                       PD_ut_mean(:), PD_vt_mean(:), PD_wt_mean(:), PD_tt_mean(:),&
                       C_uu_mean(:), C_vv_mean(:), C_ww_mean(:), &
                       C_uv_mean(:), C_uw_mean(:), C_vw_mean(:), &
                       C_ut_mean(:), C_vt_mean(:), C_wt_mean(:), C_tt_mean(:), Difv_ut_tot(:), Difv_vt_tot(:), Difv_wt_tot(:)

end module  
