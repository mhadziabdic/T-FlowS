!==============================================================================!
  subroutine Allocate_Variables(grid)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Allocate memory for velocity components ...
  call Var_Mod_Allocate_Solution("U", u, grid)
  call Var_Mod_Allocate_Solution("V", v, grid)
  call Var_Mod_Allocate_Solution("W", w, grid)

  ! ... and their gradients
  call Var_Mod_Allocate_Gradients(u)
  call Var_Mod_Allocate_Gradients(v)
  call Var_Mod_Allocate_Gradients(w)

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only("P",  p,  grid)
  call Var_Mod_Allocate_New_Only("PP", pp, grid)

  ! Pressure gradients are needed too 
  call Var_Mod_Allocate_Gradients(p)

  ! It is always calling this - probably not needed
  call Var_Mod_Allocate_Statistics(u)
  call Var_Mod_Allocate_Statistics(v)
  call Var_Mod_Allocate_Statistics(w)
  call Var_Mod_Allocate_Statistics(p)

  allocate (PHI1x(-grid % n_bnd_cells:grid % n_cells)); PHI1x=0.
  allocate (PHI1y(-grid % n_bnd_cells:grid % n_cells)); PHI1y=0.
  allocate (PHI1z(-grid % n_bnd_cells:grid % n_cells)); PHI1z=0.

  allocate (PHI2x(-grid % n_bnd_cells:grid % n_cells)); PHI2x=0.
  allocate (PHI2y(-grid % n_bnd_cells:grid % n_cells)); PHI2y=0.
  allocate (PHI2z(-grid % n_bnd_cells:grid % n_cells)); PHI2z=0.

  allocate (PHI3x(-grid % n_bnd_cells:grid % n_cells)); PHI3x=0.
  allocate (PHI3y(-grid % n_bnd_cells:grid % n_cells)); PHI3y=0.
  allocate (PHI3z(-grid % n_bnd_cells:grid % n_cells)); PHI3z=0.

  allocate (PHI4x(-grid % n_bnd_cells:grid % n_cells)); PHI4x=0.
  allocate (PHI4y(-grid % n_bnd_cells:grid % n_cells)); PHI4y=0.
  allocate (PHI4z(-grid % n_bnd_cells:grid % n_cells)); PHI4z=0.

  allocate (PHIside(grid % n_faces)); PHIside=0.

  allocate (phi_max(-grid % n_bnd_cells:grid % n_cells)); phi_max=0.
  allocate (phi_min(-grid % n_bnd_cells:grid % n_cells)); phi_min=0.

  allocate (G(6,grid % n_cells)); G=0

  allocate (Flux(grid % n_faces));     Flux=0.

  allocate (PdropX(grid % n_materials)); PdropX=0.0
  allocate (PdropY(grid % n_materials)); PdropY=0.0 
  allocate (PdropZ(grid % n_materials)); PdropZ=0.0 

  allocate (Utau(grid % n_materials));   Utau=0.0
  allocate (Vtau(grid % n_materials));   Vtau=0.0
  allocate (Wtau(grid % n_materials));   Wtau=0.0

  allocate (FLUXx(grid % n_materials));  FLUXx=0.0
  allocate (FLUXy(grid % n_materials));  FLUXy=0.0
  allocate (FLUXz(grid % n_materials));  FLUXz=0.0

  allocate (FLUXoX(grid % n_materials)); FLUXoX=0.0
  allocate (FLUXoY(grid % n_materials)); FLUXoY=0.0
  allocate (FLUXoZ(grid % n_materials)); FLUXoZ=0.0

  allocate (Ubulk(grid % n_materials));  Ubulk=0.0
  allocate (Vbulk(grid % n_materials));  Vbulk=0.0
  allocate (Wbulk(grid % n_materials));  Wbulk=0.0

  allocate (MassIn(grid % n_materials)); MassIn=0.0
  allocate (MasOut(grid % n_materials)); MasOut=0.0

  allocate (BadForG(grid % n_cells));  BadForG = .false.
  allocate (NumGood(grid % n_cells));  NumGood = 0          
  allocate (NumNeig(grid % n_cells));  NumNeig = 0         

  allocate (near(-grid % n_bnd_cells:grid % n_cells));  near  = 0.
  allocate (VISwall(-grid % n_bnd_cells:grid % n_cells)); VISwall =0.0

  ! For solution of temperature
  if(HOT==YES) then
    call Var_Mod_Allocate_Solution("T", t, grid)
    allocate (CONwall(-grid % n_bnd_cells:grid % n_cells)); CONwall =0.0
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(SIMULA==EBM.or.SIMULA==HJ) then
    allocate (VAR1x(-grid % n_bnd_cells:grid % n_cells)); VAR1x=0.
    allocate (VAR1y(-grid % n_bnd_cells:grid % n_cells)); VAR1y=0.
    allocate (VAR1z(-grid % n_bnd_cells:grid % n_cells)); VAR1z=0.

    allocate (VAR2x(-grid % n_bnd_cells:grid % n_cells)); VAR2x=0.
    allocate (VAR2y(-grid % n_bnd_cells:grid % n_cells)); VAR2y=0.
    allocate (VAR2z(-grid % n_bnd_cells:grid % n_cells)); VAR2z=0.
    if(URANS == YES) then    
      allocate (VAR10x(-grid % n_bnd_cells:grid % n_cells)); VAR10x=0.
      allocate (VAR10y(-grid % n_bnd_cells:grid % n_cells)); VAR10y=0.
      allocate (VAR10z(-grid % n_bnd_cells:grid % n_cells)); VAR10z=0.

      allocate (VAR11x(-grid % n_bnd_cells:grid % n_cells)); VAR11x=0.
      allocate (VAR11y(-grid % n_bnd_cells:grid % n_cells)); VAR11y=0.
      allocate (VAR11z(-grid % n_bnd_cells:grid % n_cells)); VAR11z=0.
    end if
    if(SIMULA == HJ) then
!      allocate (VAR3x(-grid % n_bnd_cells:grid % n_cells)); VAR3x=0.
!      allocate (VAR3y(-grid % n_bnd_cells:grid % n_cells)); VAR3y=0.
!      allocate (VAR3z(-grid % n_bnd_cells:grid % n_cells)); VAR3z=0.

!      allocate (VAR4x(-grid % n_bnd_cells:grid % n_cells)); VAR4x=0.
!      allocate (VAR4y(-grid % n_bnd_cells:grid % n_cells)); VAR4y=0.
!      allocate (VAR4z(-grid % n_bnd_cells:grid % n_cells)); VAR4z=0.

!      allocate (VAR5x(-grid % n_bnd_cells:grid % n_cells)); VAR5x=0.
!      allocate (VAR5y(-grid % n_bnd_cells:grid % n_cells)); VAR5y=0.
!      allocate (VAR5z(-grid % n_bnd_cells:grid % n_cells)); VAR5z=0.

      allocate (VAR6x(-grid % n_bnd_cells:grid % n_cells)); VAR6x=0.
      allocate (VAR6y(-grid % n_bnd_cells:grid % n_cells)); VAR6y=0.
      allocate (VAR6z(-grid % n_bnd_cells:grid % n_cells)); VAR6z=0.

      allocate (VAR7x(-grid % n_bnd_cells:grid % n_cells)); VAR7x=0.
      allocate (VAR7y(-grid % n_bnd_cells:grid % n_cells)); VAR7y=0.
      allocate (VAR7z(-grid % n_bnd_cells:grid % n_cells)); VAR7z=0.

      allocate (Eps_tot(-grid % n_bnd_cells:grid % n_cells)); Eps_tot=0.

!      allocate (VAR8x(-grid % n_bnd_cells:grid % n_cells)); VAR8x=0.
!      allocate (VAR8y(-grid % n_bnd_cells:grid % n_cells)); VAR8y=0.
!      allocate (VAR8z(-grid % n_bnd_cells:grid % n_cells)); VAR8z=0.

!      allocate (VAR9x(-grid % n_bnd_cells:grid % n_cells)); VAR9x=0.
!      allocate (VAR9y(-grid % n_bnd_cells:grid % n_cells)); VAR9y=0.
!      allocate (VAR9z(-grid % n_bnd_cells:grid % n_cells)); VAR9z=0.

!      allocate (VAR12x(-grid % n_bnd_cells:grid % n_cells)); VAR12x=0.
!      allocate (VAR12y(-grid % n_bnd_cells:grid % n_cells)); VAR12y=0.
!      allocate (VAR12z(-grid % n_bnd_cells:grid % n_cells)); VAR12z=0.
    end if  ! SIMULA == HJ

    ! Reynolds stresses
    call Var_Mod_Allocate_Solution("UU", uu, grid)
    call Var_Mod_Allocate_Solution("VV", vv, grid)
    call Var_Mod_Allocate_Solution("WW", ww, grid)
    call Var_Mod_Allocate_Solution("UV", uv, grid)
    call Var_Mod_Allocate_Solution("UW", uw, grid)
    call Var_Mod_Allocate_Solution("VW", vw, grid)

    call Var_Mod_Allocate_New_Only("KIN", kin, grid)
    call Var_Mod_Allocate_Solution("EPS", eps, grid)

    ! Time scale, length scale and production
    allocate (Tsc(-grid % n_bnd_cells:grid % n_cells));     Tsc = 0.0
    allocate (Lsc(-grid % n_bnd_cells:grid % n_cells));     Lsc = 0.0
    allocate (Pk(-grid % n_bnd_cells:grid % n_cells));      Pk  = 0.0

    if(SIMULA==EBM) then
      call Var_Mod_Allocate_Solution("F22", f22, grid)
    else
      call Var_Mod_Allocate_New_Only("F22", f22, grid)
    end if

    if(URANS == YES) then
      call Var_Mod_Allocate_Statistics(uu)
      call Var_Mod_Allocate_Statistics(vv)
      call Var_Mod_Allocate_Statistics(ww)
      call Var_Mod_Allocate_Statistics(uv)
      call Var_Mod_Allocate_Statistics(uw)
      call Var_Mod_Allocate_Statistics(vw)
      call Var_Mod_Allocate_Statistics(kin)
    end if
  end if  ! SIMULA == EBM or HJ

  ! Variables for Rans models
  if(SIMULA==K_EPS.or.SIMULA == HYB_PITM) then

    call Var_Mod_Allocate_Solution("KIN", kin, grid)
    call Var_Mod_Allocate_Solution("EPS", eps, grid)

    allocate (Uf(-grid % n_bnd_cells:grid % n_cells));      Uf    =0.0
    allocate (Ufmean(-grid % n_bnd_cells:grid % n_cells));  Ufmean=0.0
    allocate (Pk(-grid % n_bnd_cells:grid % n_cells));      Pk    =0.0
    allocate (Ynd(-grid % n_bnd_cells:grid % n_cells));     Ynd   =0.0

    if(URANS == YES) then
      allocate (Kin % mean(grid % n_cells));   Kin % mean=0.
      allocate (Eps % mean(grid % n_cells));   Eps % mean=0.
    end if
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
    call Var_Mod_Allocate_Solution("KIN", kin, grid)
    call Var_Mod_Allocate_Solution("EPS", eps, grid)
    call Var_Mod_Allocate_Solution("V^2", v_2, grid)
    call Var_Mod_Allocate_Solution("F22", f22, grid)

    allocate (Tsc(-grid % n_bnd_cells:grid % n_cells));     Tsc   =0.0
    allocate (Lsc(-grid % n_bnd_cells:grid % n_cells));     Lsc   =0.0
    allocate (Uf(-grid % n_bnd_cells:grid % n_cells));      Uf    =0.0
    allocate (Ufmean(-grid % n_bnd_cells:grid % n_cells));  Ufmean=0.0  
    allocate (Pk(-grid % n_bnd_cells:grid % n_cells));      Pk    =0.0  
    allocate (Ynd(-grid % n_bnd_cells:grid % n_cells));     Ynd   =0.0

    if(URANS == YES) then
      allocate (Kin % mean(grid % n_cells));   Kin % mean=0.
      allocate (Eps % mean(grid % n_cells));   Eps % mean=0.
      allocate (f22 % mean(grid % n_cells));   f22 % mean=0.
      allocate (v_2 % mean(grid % n_cells));   v_2 % mean=0.
    end if

    if(BUOY==YES) then
      call Var_Mod_Allocate_Solution("TT", tt, grid)
      call Var_Mod_Allocate_Statistics(tt)
      allocate (Gbuoy(-grid % n_bnd_cells:grid % n_cells));  Gbuoy=0. ! XXXXX 6 Jun 2014
      allocate (buoyBeta(-grid % n_bnd_cells:grid % n_cells));  buoyBeta=0.  ! XXXXX 5 Jul 2014
      allocate (Pbuoy(-grid % n_bnd_cells:grid % n_cells));  Pbuoy=0.  ! XXXXX 5 Jul 2014
      allocate (Kin%mean(-grid % n_bnd_cells:grid % n_cells));  Kin%mean=0.  ! XXXXX 5 Jul 2014
      allocate (Eps%mean(-grid % n_bnd_cells:grid % n_cells));  Eps%mean=0.  ! XXXXX 5 Jul 2014
      allocate (Ptt(-grid % n_bnd_cells:grid % n_cells));    Ptt=0.
    end if
  end if                    

  if(SIMULA == HYB_ZETA) then
    allocate (VISt_sgs(-grid % n_bnd_cells:grid % n_cells));  VISt_sgs=0.
    allocate (VISt_eff(-grid % n_bnd_cells:grid % n_cells));  VISt_eff=0.
  end if

  if(SIMULA == DES_SPA) then
    allocate (Ksgs(-grid % n_bnd_cells:grid % n_cells));  Ksgs=0.
  end if

  if(SIMULA == SPA_ALL.or.SIMULA == DES_SPA) then
    call Var_Mod_Allocate_Solution("VIS", vis, grid)
  end if

  if(SIMULA == DES_SPA) then
    allocate (VIS % mean(grid % n_cells));   VIS % mean=0.
  end if

  ! Variables defined in les_mod.h90:
  if(SIMULA == LES.or.SIMULA==HYB_ZETA) then
    if(MODE == WALE) then 
      allocate (WALEv(-grid % n_bnd_cells:grid % n_cells));  WALEv =0.
    end if
    if(MODE == DYN) then 
      allocate (U % filt(-grid % n_bnd_cells:grid % n_cells));  U % filt =0.
      allocate (V % filt(-grid % n_bnd_cells:grid % n_cells));  V % filt =0.
      allocate (W % filt(-grid % n_bnd_cells:grid % n_cells));  W % filt =0.
   
      allocate (Cdyn(-grid % n_bnd_cells:grid % n_cells)); Cdyn = 0
      allocate(UUf(grid % n_cells));   UUf = 0.0
      allocate(VVf(grid % n_cells));   VVf = 0.0
      allocate(WWf(grid % n_cells));   WWf = 0.0
      allocate(UVf(grid % n_cells));   UVf = 0.0
      allocate(UWf(grid % n_cells));   UWf = 0.0
      allocate(VWf(grid % n_cells));   VWf = 0.0

      allocate(M11f(grid % n_cells));   M11f = 0.0
      allocate(M22f(grid % n_cells));   M22f = 0.0
      allocate(M33f(grid % n_cells));   M33f = 0.0
      allocate(M12f(grid % n_cells));   M12f = 0.0
      allocate(M13f(grid % n_cells));   M13f = 0.0
      allocate(M23f(grid % n_cells));   M23f = 0.0
    end if 
    allocate(ShearTest(-grid % n_bnd_cells:grid % n_cells));   ShearTest = 0.0
    allocate (Ksgs(-grid % n_bnd_cells:grid % n_cells));  Ksgs=0.
    allocate (Cdyn_mean(-grid % n_bnd_cells:grid % n_cells)); Cdyn_mean = 0
  end if

  if(SIMULA == LES.or.SIMULA==DNS.or.SIMULA==DES_SPA) then
    allocate (uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean=0.
    allocate (vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean=0.
    allocate (ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean=0.
    allocate (uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean=0.
    allocate (uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean=0.
    allocate (vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean=0.

    allocate (uuu % mean(-grid % n_bnd_cells:grid % n_cells)); uuu % mean=0.
    allocate (uuv % mean(-grid % n_bnd_cells:grid % n_cells)); uuv % mean=0.
    allocate (uuw % mean(-grid % n_bnd_cells:grid % n_cells)); uuw % mean=0.

    allocate (vvu % mean(-grid % n_bnd_cells:grid % n_cells)); vvu % mean=0.
    allocate (vvv % mean(-grid % n_bnd_cells:grid % n_cells)); vvv % mean=0.
    allocate (vvw % mean(-grid % n_bnd_cells:grid % n_cells)); vvw % mean=0.

    allocate (wwu % mean(-grid % n_bnd_cells:grid % n_cells)); wwu % mean=0.
    allocate (wwv % mean(-grid % n_bnd_cells:grid % n_cells)); wwv % mean=0.
    allocate (www % mean(-grid % n_bnd_cells:grid % n_cells)); www % mean=0.

    allocate (uwu % mean(-grid % n_bnd_cells:grid % n_cells)); uwu % mean=0.
    allocate (uwv % mean(-grid % n_bnd_cells:grid % n_cells)); uwv % mean=0.
    allocate (uww % mean(-grid % n_bnd_cells:grid % n_cells)); uww % mean=0.

    if(BUDG==YES) then
      allocate (uu % n(-grid % n_bnd_cells:grid % n_cells)); uu % n=0.
      allocate (vv % n(-grid % n_bnd_cells:grid % n_cells)); vv % n=0.
      allocate (ww % n(-grid % n_bnd_cells:grid % n_cells)); ww % n=0.
      allocate (uv % n(-grid % n_bnd_cells:grid % n_cells)); uv % n=0.
      allocate (uw % n(-grid % n_bnd_cells:grid % n_cells)); uw % n=0.
      allocate (vw % n(-grid % n_bnd_cells:grid % n_cells)); vw % n=0.

      allocate (U % fluc(-grid % n_bnd_cells:grid % n_cells));  U % fluc =0.
      allocate (V % fluc(-grid % n_bnd_cells:grid % n_cells));  V % fluc =0.
      allocate (W % fluc(-grid % n_bnd_cells:grid % n_cells));  W % fluc =0.
      allocate (P % fluc(-grid % n_bnd_cells:grid % n_cells));  P % fluc =0.

      allocate (Puu_mean(1:grid % n_cells));  Puu_mean =0.
      allocate (Pvv_mean(1:grid % n_cells));  Pvv_mean =0.
      allocate (Pww_mean(1:grid % n_cells));  Pww_mean =0.
      allocate (Puv_mean(1:grid % n_cells));  Puv_mean =0.
      allocate (Puw_mean(1:grid % n_cells));  Puw_mean =0.
      allocate (Pvw_mean(1:grid % n_cells));  Pvw_mean =0.

      allocate (Diss_uu_mean(1:grid % n_cells)); Diss_uu_mean  =0.
      allocate (Diss_vv_mean(1:grid % n_cells)); Diss_vv_mean  =0.
      allocate (Diss_ww_mean(1:grid % n_cells)); Diss_ww_mean  =0.
      allocate (Diss_uv_mean(1:grid % n_cells)); Diss_uv_mean  =0.
      allocate (Diss_uw_mean(1:grid % n_cells)); Diss_uw_mean  =0.
      allocate (Diss_vw_mean(1:grid % n_cells)); Diss_vw_mean  =0.

      allocate (Diss_sgs_mean(1:grid % n_cells)); Diss_sgs_mean  =0.

      allocate (PHI5x(-grid % n_bnd_cells:grid % n_cells)); PHI5x=0.
      allocate (PHI5y(-grid % n_bnd_cells:grid % n_cells)); PHI5y=0.
      allocate (PHI5z(-grid % n_bnd_cells:grid % n_cells)); PHI5z=0.

      allocate (PHI6x(-grid % n_bnd_cells:grid % n_cells)); PHI6x=0.
      allocate (PHI6y(-grid % n_bnd_cells:grid % n_cells)); PHI6y=0.
      allocate (PHI6z(-grid % n_bnd_cells:grid % n_cells)); PHI6z=0.

      allocate (PHI7x(-grid % n_bnd_cells:grid % n_cells)); PHI7x=0.
      allocate (PHI7y(-grid % n_bnd_cells:grid % n_cells)); PHI7y=0.
      allocate (PHI7z(-grid % n_bnd_cells:grid % n_cells)); PHI7z=0.

      allocate (PHI8x(-grid % n_bnd_cells:grid % n_cells)); PHI8x=0.
      allocate (PHI8y(-grid % n_bnd_cells:grid % n_cells)); PHI8y=0.
      allocate (PHI8z(-grid % n_bnd_cells:grid % n_cells)); PHI8z=0.

      allocate (PHI9x(-grid % n_bnd_cells:grid % n_cells)); PHI9x=0.
      allocate (PHI9y(-grid % n_bnd_cells:grid % n_cells)); PHI9y=0.
      allocate (PHI9z(-grid % n_bnd_cells:grid % n_cells)); PHI9z=0.

      allocate (PHI10x(-grid % n_bnd_cells:grid % n_cells)); PHI10x=0.
      allocate (PHI10y(-grid % n_bnd_cells:grid % n_cells)); PHI10y=0.
      allocate (PHI10z(-grid % n_bnd_cells:grid % n_cells)); PHI10z=0.

      if(HOT==YES) then
        allocate (Put_mean(1:grid % n_cells));  Put_mean =0.
        allocate (Pvt_mean(1:grid % n_cells));  Pvt_mean =0.
        allocate (Pwt_mean(1:grid % n_cells));  Pwt_mean =0.
        allocate (Ptt_mean(1:grid % n_cells));  Ptt_mean =0.
        allocate (Difv_ut_tot(1:grid % n_cells));   Difv_ut_tot=0.
        allocate (Difv_vt_tot(1:grid % n_cells));   Difv_vt_tot=0.
        allocate (Difv_wt_tot(1:grid % n_cells));   Difv_wt_tot=0.
        allocate (Diss_ut_mean(1:grid % n_cells)); Diss_ut_mean  =0.
        allocate (Diss_vt_mean(1:grid % n_cells)); Diss_vt_mean  =0.
        allocate (Diss_wt_mean(1:grid % n_cells)); Diss_wt_mean  =0.
        allocate (Diss_tt_mean(1:grid % n_cells)); Diss_tt_mean  =0.
        allocate (Dift_ut_mean(1:grid % n_cells));   Dift_ut_mean=0.
        allocate (Dift_vt_mean(1:grid % n_cells));   Dift_vt_mean=0.
        allocate (Dift_wt_mean(1:grid % n_cells));   Dift_wt_mean=0.
        allocate (Dift_tt_mean(1:grid % n_cells));   Dift_tt_mean=0.
        allocate (Difv_ut_mean(1:grid % n_cells));   Difv_ut_mean=0.
        allocate (Difv_vt_mean(1:grid % n_cells));   Difv_vt_mean=0.
        allocate (Difv_wt_mean(1:grid % n_cells));   Difv_wt_mean=0.
        allocate (Difv_tt_mean(1:grid % n_cells));   Difv_tt_mean=0.
        allocate (C_ut_mean(1:grid % n_cells));   C_ut_mean=0.
        allocate (C_vt_mean(1:grid % n_cells));   C_vt_mean=0.
        allocate (C_wt_mean(1:grid % n_cells));   C_wt_mean=0.
        allocate (C_tt_mean(1:grid % n_cells));   C_tt_mean=0.
        allocate (PD_ut_mean(-grid % n_bnd_cells:grid % n_cells));   PD_ut_mean=0.
        allocate (PD_vt_mean(-grid % n_bnd_cells:grid % n_cells));   PD_vt_mean=0.
        allocate (PD_wt_mean(-grid % n_bnd_cells:grid % n_cells));   PD_wt_mean=0.
        allocate (PR_ut_mean(-grid % n_bnd_cells:grid % n_cells));   PR_ut_mean=0.
        allocate (PR_vt_mean(-grid % n_bnd_cells:grid % n_cells));   PR_vt_mean=0.
        allocate (PR_wt_mean(-grid % n_bnd_cells:grid % n_cells));   PR_wt_mean=0.
        allocate (T % fluc(-grid % n_bnd_cells:grid % n_cells));  T % fluc=0.
      end if

      allocate (Dift_uu_mean(1:grid % n_cells));   Dift_uu_mean=0.
      allocate (Dift_vv_mean(1:grid % n_cells));   Dift_vv_mean=0.
      allocate (Dift_ww_mean(1:grid % n_cells));   Dift_ww_mean=0.
      allocate (Dift_uv_mean(1:grid % n_cells));   Dift_uv_mean=0.
      allocate (Dift_uw_mean(1:grid % n_cells));   Dift_uw_mean=0.
      allocate (Dift_vw_mean(1:grid % n_cells));   Dift_vw_mean=0.

      allocate (Difv_uu_mean(1:grid % n_cells));   Difv_uu_mean=0.
      allocate (Difv_vv_mean(1:grid % n_cells));   Difv_vv_mean=0.
      allocate (Difv_ww_mean(1:grid % n_cells));   Difv_ww_mean=0.
      allocate (Difv_uv_mean(1:grid % n_cells));   Difv_uv_mean=0.
      allocate (Difv_uw_mean(1:grid % n_cells));   Difv_uw_mean=0.
      allocate (Difv_vw_mean(1:grid % n_cells));   Difv_vw_mean=0.

      allocate (C_uu_mean(1:grid % n_cells));   C_uu_mean=0.
      allocate (C_vv_mean(1:grid % n_cells));   C_vv_mean=0.
      allocate (C_ww_mean(1:grid % n_cells));   C_ww_mean=0.
      allocate (C_uv_mean(1:grid % n_cells));   C_uv_mean=0.
      allocate (C_uw_mean(1:grid % n_cells));   C_uw_mean=0.
      allocate (C_vw_mean(1:grid % n_cells));   C_vw_mean=0.

      allocate (PD_uu_mean(-grid % n_bnd_cells:grid % n_cells));   PD_uu_mean=0.
      allocate (PD_vv_mean(-grid % n_bnd_cells:grid % n_cells));   PD_vv_mean=0.
      allocate (PD_ww_mean(-grid % n_bnd_cells:grid % n_cells));   PD_ww_mean=0.
      allocate (PD_uv_mean(-grid % n_bnd_cells:grid % n_cells));   PD_uv_mean=0.
      allocate (PD_uw_mean(-grid % n_bnd_cells:grid % n_cells));   PD_uw_mean=0.
      allocate (PD_vw_mean(-grid % n_bnd_cells:grid % n_cells));   PD_vw_mean=0.

      allocate (PR_uu_mean(-grid % n_bnd_cells:grid % n_cells));   PR_uu_mean=0.
      allocate (PR_vv_mean(-grid % n_bnd_cells:grid % n_cells));   PR_vv_mean=0.
      allocate (PR_ww_mean(-grid % n_bnd_cells:grid % n_cells));   PR_ww_mean=0.
      allocate (PR_uv_mean(-grid % n_bnd_cells:grid % n_cells));   PR_uv_mean=0.
      allocate (PR_uw_mean(-grid % n_bnd_cells:grid % n_cells));   PR_uw_mean=0.
      allocate (PR_vw_mean(-grid % n_bnd_cells:grid % n_cells));   PR_vw_mean=0.
    end if

    allocate(VISt_mean(grid % n_cells)); VISt_mean = 0.0
    allocate (ShearMean(grid % n_cells));  ShearMean=0.

    if(HOT==YES) then
      allocate (T % mean(-grid % n_bnd_cells:grid % n_cells));  T % mean=0.
      allocate (TT % mean(-grid % n_bnd_cells:grid % n_cells)); TT % mean=0.
      allocate (uT % mean(-grid % n_bnd_cells:grid % n_cells)); uT % mean=0.
      allocate (vT % mean(-grid % n_bnd_cells:grid % n_cells)); vT % mean=0.
      allocate (wT % mean(-grid % n_bnd_cells:grid % n_cells)); wT % mean=0.
    end if
  end if

  if(SIMULA == HYB_ZETA.or.SIMULA==HYB_PITM) then
    allocate (uu % mean(-grid % n_bnd_cells:grid % n_cells)); uu % mean=0.
    allocate (vv % mean(-grid % n_bnd_cells:grid % n_cells)); vv % mean=0.
    allocate (ww % mean(-grid % n_bnd_cells:grid % n_cells)); ww % mean=0.
    allocate (uv % mean(-grid % n_bnd_cells:grid % n_cells)); uv % mean=0.
    allocate (uw % mean(-grid % n_bnd_cells:grid % n_cells)); uw % mean=0.
    allocate (vw % mean(-grid % n_bnd_cells:grid % n_cells)); vw % mean=0.

    allocate(VISt_mean(grid % n_cells)); VISt_mean = 0.0
    allocate (ShearMean(grid % n_cells));  ShearMean=0.
    if(HOT==YES) then
      allocate (T % mean(-grid % n_bnd_cells:grid % n_cells));  T % mean=0.
      allocate (TT % mean(-grid % n_bnd_cells:grid % n_cells)); TT % mean=0.
      allocate (uT % mean(-grid % n_bnd_cells:grid % n_cells)); uT % mean=0.
      allocate (vT % mean(-grid % n_bnd_cells:grid % n_cells)); vT % mean=0.
      allocate (wT % mean(-grid % n_bnd_cells:grid % n_cells)); wT % mean=0.
    end if
  end if
 
  allocate (VISt(-grid % n_bnd_cells:grid % n_cells)); VISt=0
  allocate (IsNearWall(grid % n_cells)); IsNearWall = .FALSE.

  allocate (Vort(-grid % n_bnd_cells:grid % n_cells));  Vort=0.
  allocate (Shear(-grid % n_bnd_cells:grid % n_cells)); Shear=0.
  allocate (TauWall(grid % n_cells));    TauWall=0.

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something !  

  end subroutine
