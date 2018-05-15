!==============================================================================!
  subroutine Allocate_Variables(grid)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
!==============================================================================!

  ! Allocate memory for velocity components ...
  call Var_Mod_Allocate_Solution('U', u, grid)
  call Var_Mod_Allocate_Solution('V', v, grid)
  call Var_Mod_Allocate_Solution('W', w, grid)

  ! ... and their gradients
  call Var_Mod_Allocate_Gradients(u)
  call Var_Mod_Allocate_Gradients(v)
  call Var_Mod_Allocate_Gradients(w)

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only('P',  p,  grid)
  call Var_Mod_Allocate_New_Only('PP', pp, grid)

  ! Pressure gradients are needed too
  call Var_Mod_Allocate_Gradients(p)

  ! It is always calling this - probably not needed
  call Var_Mod_Allocate_Statistics(u)
  call Var_Mod_Allocate_Statistics(v)
  call Var_Mod_Allocate_Statistics(w)
  call Var_Mod_Allocate_Statistics(p)

  allocate(phi_face(grid % n_faces)); phi_face = 0.

  allocate(phi_max(-grid % n_bnd_cells:grid % n_cells)); phi_max = 0.
  allocate(phi_min(-grid % n_bnd_cells:grid % n_cells)); phi_min = 0.

  allocate(flux(grid % n_faces));     Flux = 0.

  allocate(Utau(grid % n_materials));   Utau = 0.
  allocate(Vtau(grid % n_materials));   Vtau = 0.
  allocate(Wtau(grid % n_materials));   Wtau = 0.

  call Grad_Mod_Allocate_Memory(grid)

  allocate(nearest_wall_cell(-grid % n_bnd_cells:grid % n_cells))
  nearest_wall_cell = 0

  allocate(vis_wall(-grid % n_bnd_cells:grid % n_cells)); vis_wall = 0.

  ! For solution of temperature
  call Control_Mod_Heat_Transfer(verbose = .true.)
  if(heat_transfer == YES) then
    call Var_Mod_Allocate_Solution('T',  t,  grid)
    call Var_Mod_Allocate_Solution('UT', ut, grid)
    call Var_Mod_Allocate_Solution('VT', vt, grid)
    call Var_Mod_Allocate_Solution('WT', wt, grid)
    allocate(con_wall(-grid % n_bnd_cells:grid % n_cells)); con_wall = 0.
  end if

  call Control_Mod_Buoyancy(verbose = .true.)

  call Control_Mod_Turbulence_Model(verbose = .true.)
  call Control_Mod_Turbulence_Model_Variant(verbose = .true.)

  ! Allocate wall distance for all models
  if(turbulence_model .ne. NONE) then
    allocate(y_plus(-grid % n_bnd_cells:grid % n_cells));  y_plus = 0.
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant == URANS) then
!
!     Should something be done here?
!
    end if
    if(turbulence_model == HANJALIC_JAKIRLIC) then
      allocate(Eps_tot(-grid % n_bnd_cells:grid % n_cells)); Eps_tot = 0.
    end if

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution('UU', uu, grid)
    call Var_Mod_Allocate_Solution('VV', vv, grid)
    call Var_Mod_Allocate_Solution('WW', ww, grid)
    call Var_Mod_Allocate_Solution('UV', uv, grid)
    call Var_Mod_Allocate_Solution('UW', uw, grid)
    call Var_Mod_Allocate_Solution('VW', vw, grid)

    call Var_Mod_Allocate_New_Only('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    ! Time scale, length scale and production
    allocate(t_scale(-grid % n_bnd_cells:grid % n_cells));  t_scale = 0.
    allocate(l_scale(-grid % n_bnd_cells:grid % n_cells));  l_scale = 0.
    allocate(p_kin  (-grid % n_bnd_cells:grid % n_cells));  p_kin   = 0.

    if(turbulence_model == REYNOLDS_STRESS_MODEL) then
      call Var_Mod_Allocate_Solution('F22', f22, grid)
      call Var_Mod_Allocate_Gradients(f22)
    else
      call Var_Mod_Allocate_New_Only('F22', f22, grid)
    end if

    if(turbulence_model_variant == URANS) then
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)
      call Var_Mod_Allocate_Statistics(kin)
    end if
  end if  ! turbulence_model == 'EBM' or 'HJ'

  ! Variables for Rans models
  if(turbulence_model == K_EPS) then
    call Var_Mod_Allocate_Solution('KIN', kin, grid)
    call Var_Mod_Allocate_Solution('EPS', eps, grid)

    allocate(u_tau     (-grid % n_bnd_cells:grid % n_cells));  u_tau      = 0.
    allocate(u_tau_mean(-grid % n_bnd_cells:grid % n_cells));  u_tau_mean = 0.
    allocate(p_kin     (-grid % n_bnd_cells:grid % n_cells));  p_kin      = 0.

    if(turbulence_model_variant == URANS) then
      allocate(kin % mean(grid % n_cells));   kin % mean = 0.
      allocate(eps % mean(grid % n_cells));   eps % mean = 0.
    end if
  end if

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    call Var_Mod_Allocate_Solution('KIN',  kin,  grid)
    call Var_Mod_Allocate_Solution('EPS',  eps,  grid)
    call Var_Mod_Allocate_Solution('ZETA', zeta, grid)
    call Var_Mod_Allocate_Solution('F22',  f22,  grid)

    allocate(t_scale   (-grid % n_bnd_cells:grid % n_cells));  t_scale    = 0.
    allocate(l_scale   (-grid % n_bnd_cells:grid % n_cells));  l_scale    = 0.
    allocate(u_tau     (-grid % n_bnd_cells:grid % n_cells));  u_tau      = 0.
    allocate(u_tau_mean(-grid % n_bnd_cells:grid % n_cells));  u_tau_mean = 0.
    allocate(p_kin     (-grid % n_bnd_cells:grid % n_cells));  p_kin      = 0.

    if(turbulence_model_variant == URANS) then
      allocate(kin  % mean(grid % n_cells));  kin  % mean = 0.
      allocate(eps  % mean(grid % n_cells));  eps  % mean = 0.
      allocate(zeta % mean(grid % n_cells));  zeta % mean = 0.
      allocate(f22  % mean(grid % n_cells));  f22  % mean = 0.
    end if

    if(buoyancy == YES) then
      call Var_Mod_Allocate_Solution('TT', tt, grid)
      call Var_Mod_Allocate_Statistics(tt)
      allocate(g_buoy   (-grid % n_bnd_cells:grid % n_cells));  g_buoy     = 0.
      allocate(buoy_beta(-grid % n_bnd_cells:grid % n_cells));  buoy_beta  = 0.
      allocate(p_buoy   (-grid % n_bnd_cells:grid % n_cells));  p_buoy     = 0.
      allocate(kin%mean (-grid % n_bnd_cells:grid % n_cells));  kin % mean = 0.
      allocate(eps%mean (-grid % n_bnd_cells:grid % n_cells));  eps % mean = 0.
      allocate(Ptt      (-grid % n_bnd_cells:grid % n_cells));  Ptt        = 0.
    end if
  end if

  if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    allocate(vis_t_sgs(-grid % n_bnd_cells:grid % n_cells));  vis_t_sgs = 0.
    allocate(vis_t_eff(-grid % n_bnd_cells:grid % n_cells));  vis_t_eff = 0.
  end if

  if(turbulence_model == DES_SPALART) then
    allocate(kin_sgs(-grid % n_bnd_cells:grid % n_cells));  kin_sgs= 0.
  end if

  if(turbulence_model == SPALART_ALLMARAS .or.  &
     turbulence_model == DES_SPALART) then
    call Var_Mod_Allocate_Solution('VIS', vis, grid)
  end if

  if(turbulence_model == DES_SPALART) then
    allocate(VIS % mean(grid % n_cells));   VIS % mean= 0.
  end if

  ! Variables defined in les_mod.h90:
  if(turbulence_model == LES .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    if(turbulence_model_variant == WALE) then
      allocate(wale_v(-grid % n_bnd_cells:grid % n_cells));  wale_v = 0.
    end if
    if(turbulence_model_variant == DYNAMIC) then
      allocate(u % filt(-grid % n_bnd_cells:grid % n_cells));  u % filt = 0.
      allocate(v % filt(-grid % n_bnd_cells:grid % n_cells));  v % filt = 0.
      allocate(w % filt(-grid % n_bnd_cells:grid % n_cells));  w % filt = 0.

      allocate(c_dyn(-grid % n_bnd_cells:grid % n_cells)); c_dyn = 0.
      allocate(UUf(grid % n_cells));   UUf = 0.
      allocate(VVf(grid % n_cells));   VVf = 0.
      allocate(WWf(grid % n_cells));   WWf = 0.
      allocate(UVf(grid % n_cells));   UVf = 0.
      allocate(UWf(grid % n_cells));   UWf = 0.
      allocate(VWf(grid % n_cells));   VWf = 0.

      allocate(M11f(grid % n_cells));   M11f = 0.
      allocate(M22f(grid % n_cells));   M22f = 0.
      allocate(M33f(grid % n_cells));   M33f = 0.
      allocate(M12f(grid % n_cells));   M12f = 0.
      allocate(M13f(grid % n_cells));   M13f = 0.
      allocate(M23f(grid % n_cells));   M23f = 0.
    end if
    allocate(shear_test(-grid % n_bnd_cells:grid % n_cells)); shear_test = 0.
    allocate(kin_sgs   (-grid % n_bnd_cells:grid % n_cells)); kin_sgs= 0.
    allocate(c_dyn_mean(-grid % n_bnd_cells:grid % n_cells)); c_dyn_mean = 0.
  end if

  if(turbulence_model == LES .or.  &
     turbulence_model == DNS .or.  &
     turbulence_model == DES_SPALART) then
    allocate(uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean= 0.
    allocate(vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean= 0.
    allocate(ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean= 0.
    allocate(uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean= 0.
    allocate(uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean= 0.
    allocate(vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean= 0.

    allocate(uuu % mean(-grid % n_bnd_cells:grid % n_cells)); uuu % mean= 0.
    allocate(uuv % mean(-grid % n_bnd_cells:grid % n_cells)); uuv % mean= 0.
    allocate(uuw % mean(-grid % n_bnd_cells:grid % n_cells)); uuw % mean= 0.

    allocate(vvu % mean(-grid % n_bnd_cells:grid % n_cells)); vvu % mean= 0.
    allocate(vvv % mean(-grid % n_bnd_cells:grid % n_cells)); vvv % mean= 0.
    allocate(vvw % mean(-grid % n_bnd_cells:grid % n_cells)); vvw % mean= 0.

    allocate(wwu % mean(-grid % n_bnd_cells:grid % n_cells)); wwu % mean= 0.
    allocate(wwv % mean(-grid % n_bnd_cells:grid % n_cells)); wwv % mean= 0.
    allocate(www % mean(-grid % n_bnd_cells:grid % n_cells)); www % mean= 0.

    allocate(uwu % mean(-grid % n_bnd_cells:grid % n_cells)); uwu % mean= 0.
    allocate(uwv % mean(-grid % n_bnd_cells:grid % n_cells)); uwv % mean= 0.
    allocate(uww % mean(-grid % n_bnd_cells:grid % n_cells)); uww % mean= 0.

    if(BUDG==YES) then
      allocate(uu % n(-grid % n_bnd_cells:grid % n_cells)); uu % n= 0.
      allocate(vv % n(-grid % n_bnd_cells:grid % n_cells)); vv % n= 0.
      allocate(ww % n(-grid % n_bnd_cells:grid % n_cells)); ww % n= 0.
      allocate(uv % n(-grid % n_bnd_cells:grid % n_cells)); uv % n= 0.
      allocate(uw % n(-grid % n_bnd_cells:grid % n_cells)); uw % n= 0.
      allocate(vw % n(-grid % n_bnd_cells:grid % n_cells)); vw % n= 0.

      allocate(u % fluc(-grid % n_bnd_cells:grid % n_cells));  u % fluc = 0.
      allocate(v % fluc(-grid % n_bnd_cells:grid % n_cells));  v % fluc = 0.
      allocate(w % fluc(-grid % n_bnd_cells:grid % n_cells));  w % fluc = 0.
      allocate(p % fluc(-grid % n_bnd_cells:grid % n_cells));  p % fluc = 0.

      allocate(Puu_mean(1:grid % n_cells));  Puu_mean = 0.
      allocate(Pvv_mean(1:grid % n_cells));  Pvv_mean = 0.
      allocate(Pww_mean(1:grid % n_cells));  Pww_mean = 0.
      allocate(Puv_mean(1:grid % n_cells));  Puv_mean = 0.
      allocate(Puw_mean(1:grid % n_cells));  Puw_mean = 0.
      allocate(Pvw_mean(1:grid % n_cells));  Pvw_mean = 0.

      allocate(Diss_uu_mean(1:grid % n_cells)); Diss_uu_mean  = 0.
      allocate(Diss_vv_mean(1:grid % n_cells)); Diss_vv_mean  = 0.
      allocate(Diss_ww_mean(1:grid % n_cells)); Diss_ww_mean  = 0.
      allocate(Diss_uv_mean(1:grid % n_cells)); Diss_uv_mean  = 0.
      allocate(Diss_uw_mean(1:grid % n_cells)); Diss_uw_mean  = 0.
      allocate(Diss_vw_mean(1:grid % n_cells)); Diss_vw_mean  = 0.

      allocate(Diss_sgs_mean(1:grid % n_cells)); Diss_sgs_mean  = 0.

      if(heat_transfer == YES) then
        allocate(Put_mean(1:grid % n_cells));  Put_mean = 0.
        allocate(Pvt_mean(1:grid % n_cells));  Pvt_mean = 0.
        allocate(Pwt_mean(1:grid % n_cells));  Pwt_mean = 0.
        allocate(Ptt_mean(1:grid % n_cells));  Ptt_mean = 0.
        allocate(Difv_ut_tot(1:grid % n_cells));   Difv_ut_tot= 0.
        allocate(Difv_vt_tot(1:grid % n_cells));   Difv_vt_tot= 0.
        allocate(Difv_wt_tot(1:grid % n_cells));   Difv_wt_tot= 0.
        allocate(Diss_ut_mean(1:grid % n_cells)); Diss_ut_mean  = 0.
        allocate(Diss_vt_mean(1:grid % n_cells)); Diss_vt_mean  = 0.
        allocate(Diss_wt_mean(1:grid % n_cells)); Diss_wt_mean  = 0.
        allocate(Diss_tt_mean(1:grid % n_cells)); Diss_tt_mean  = 0.
        allocate(Dift_ut_mean(1:grid % n_cells));   Dift_ut_mean = 0.
        allocate(Dift_vt_mean(1:grid % n_cells));   Dift_vt_mean = 0.
        allocate(Dift_wt_mean(1:grid % n_cells));   Dift_wt_mean = 0.
        allocate(Dift_tt_mean(1:grid % n_cells));   Dift_tt_mean = 0.
        allocate(Difv_ut_mean(1:grid % n_cells));   Difv_ut_mean = 0.
        allocate(Difv_vt_mean(1:grid % n_cells));   Difv_vt_mean = 0.
        allocate(Difv_wt_mean(1:grid % n_cells));   Difv_wt_mean = 0.
        allocate(Difv_tt_mean(1:grid % n_cells));   Difv_tt_mean = 0.
        allocate(C_ut_mean(1:grid % n_cells));   C_ut_mean = 0.
        allocate(C_vt_mean(1:grid % n_cells));   C_vt_mean = 0.
        allocate(C_wt_mean(1:grid % n_cells));   C_wt_mean = 0.
        allocate(C_tt_mean(1:grid % n_cells));   C_tt_mean = 0.
        allocate(PD_ut_mean(-grid % n_bnd_cells:grid % n_cells)); PD_ut_mean= 0.
        allocate(PD_vt_mean(-grid % n_bnd_cells:grid % n_cells)); PD_vt_mean= 0.
        allocate(PD_wt_mean(-grid % n_bnd_cells:grid % n_cells)); PD_wt_mean= 0.
        allocate(PR_ut_mean(-grid % n_bnd_cells:grid % n_cells)); PR_ut_mean= 0.
        allocate(PR_vt_mean(-grid % n_bnd_cells:grid % n_cells)); PR_vt_mean= 0.
        allocate(PR_wt_mean(-grid % n_bnd_cells:grid % n_cells)); PR_wt_mean= 0.
        allocate(T % fluc(-grid % n_bnd_cells:grid % n_cells)); T % fluc= 0.
      end if

      allocate(Dift_uu_mean(1:grid % n_cells));   Dift_uu_mean = 0.
      allocate(Dift_vv_mean(1:grid % n_cells));   Dift_vv_mean = 0.
      allocate(Dift_ww_mean(1:grid % n_cells));   Dift_ww_mean = 0.
      allocate(Dift_uv_mean(1:grid % n_cells));   Dift_uv_mean = 0.
      allocate(Dift_uw_mean(1:grid % n_cells));   Dift_uw_mean = 0.
      allocate(Dift_vw_mean(1:grid % n_cells));   Dift_vw_mean = 0.

      allocate(Difv_uu_mean(1:grid % n_cells));   Difv_uu_mean = 0.
      allocate(Difv_vv_mean(1:grid % n_cells));   Difv_vv_mean = 0.
      allocate(Difv_ww_mean(1:grid % n_cells));   Difv_ww_mean = 0.
      allocate(Difv_uv_mean(1:grid % n_cells));   Difv_uv_mean = 0.
      allocate(Difv_uw_mean(1:grid % n_cells));   Difv_uw_mean = 0.
      allocate(Difv_vw_mean(1:grid % n_cells));   Difv_vw_mean = 0.

      allocate(C_uu_mean(1:grid % n_cells));   C_uu_mean = 0.
      allocate(C_vv_mean(1:grid % n_cells));   C_vv_mean = 0.
      allocate(C_ww_mean(1:grid % n_cells));   C_ww_mean = 0.
      allocate(C_uv_mean(1:grid % n_cells));   C_uv_mean = 0.
      allocate(C_uw_mean(1:grid % n_cells));   C_uw_mean = 0.
      allocate(C_vw_mean(1:grid % n_cells));   C_vw_mean = 0.

      allocate(PD_uu_mean(-grid % n_bnd_cells:grid % n_cells)); PD_uu_mean = 0.
      allocate(PD_vv_mean(-grid % n_bnd_cells:grid % n_cells)); PD_vv_mean = 0.
      allocate(PD_ww_mean(-grid % n_bnd_cells:grid % n_cells)); PD_ww_mean = 0.
      allocate(PD_uv_mean(-grid % n_bnd_cells:grid % n_cells)); PD_uv_mean = 0.
      allocate(PD_uw_mean(-grid % n_bnd_cells:grid % n_cells)); PD_uw_mean = 0.
      allocate(PD_vw_mean(-grid % n_bnd_cells:grid % n_cells)); PD_vw_mean = 0.

      allocate(PR_uu_mean(-grid % n_bnd_cells:grid % n_cells)); PR_uu_mean = 0.
      allocate(PR_vv_mean(-grid % n_bnd_cells:grid % n_cells)); PR_vv_mean = 0.
      allocate(PR_ww_mean(-grid % n_bnd_cells:grid % n_cells)); PR_ww_mean = 0.
      allocate(PR_uv_mean(-grid % n_bnd_cells:grid % n_cells)); PR_uv_mean = 0.
      allocate(PR_uw_mean(-grid % n_bnd_cells:grid % n_cells)); PR_uw_mean = 0.
      allocate(PR_vw_mean(-grid % n_bnd_cells:grid % n_cells)); PR_vw_mean = 0.
    end if

    allocate(vis_t_mean(grid % n_cells)); vis_t_mean = 0.
    allocate(shear_mean(grid % n_cells)); shear_mean = 0.

    if(heat_transfer == YES) then
      allocate(t %  mean(-grid % n_bnd_cells:grid % n_cells)); t  % mean = 0.
      allocate(tt % mean(-grid % n_bnd_cells:grid % n_cells)); tt % mean = 0.
      allocate(ut % mean(-grid % n_bnd_cells:grid % n_cells)); ut % mean = 0.
      allocate(vt % mean(-grid % n_bnd_cells:grid % n_cells)); vt % mean = 0.
      allocate(wt % mean(-grid % n_bnd_cells:grid % n_cells)); wt % mean = 0.
    end if
  end if

  if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    allocate(uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean = 0.
    allocate(vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean = 0.
    allocate(ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean = 0.
    allocate(uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean = 0.
    allocate(uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean = 0.
    allocate(vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean = 0.

    allocate(vis_t_mean(grid % n_cells));  vis_t_mean = 0.
    allocate(shear_mean(grid % n_cells));  shear_mean = 0.
    if(heat_transfer == YES) then
      allocate(t  % mean(-grid % n_bnd_cells:grid % n_cells)); t  % mean = 0.
      allocate(tt % mean(-grid % n_bnd_cells:grid % n_cells)); tt % mean = 0.
      allocate(ut % mean(-grid % n_bnd_cells:grid % n_cells)); ut % mean = 0.
      allocate(vt % mean(-grid % n_bnd_cells:grid % n_cells)); vt % mean = 0.
      allocate(wt % mean(-grid % n_bnd_cells:grid % n_cells)); wt % mean = 0.
    end if
  end if

  allocate(vis_t   (-grid % n_bnd_cells:grid % n_cells));  vis_t    = 0.
  allocate(vort    (-grid % n_bnd_cells:grid % n_cells));  vort     = 0.
  allocate(shear   (-grid % n_bnd_cells:grid % n_cells));  shear    = 0.
  allocate(tau_wall(grid % n_cells));                      tau_wall = 0.

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something !

  end subroutine
