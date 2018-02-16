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
  use Parameters_Mod
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

  ! Mean velocity and fluctuating velocities for DNS and LES
  if(SIMULA==LES) then
    call Var_Mod_Allocate_Statistics(u)
    call Var_Mod_Allocate_Statistics(v)
    call Var_Mod_Allocate_Statistics(w)
    call Var_Mod_Allocate_Statistics(p)
  end if

  if(SIMULA==HYB_ZETA.or.SIMULA==DES_SPA.or.&
     SIMULA==DNS.or.SIMULA==HYB_PITM.or.URANS==YES) then
    allocate (U % mean(grid % n_cells));   U % mean=0.
    allocate (V % mean(grid % n_cells));   V % mean=0.
    allocate (W % mean(grid % n_cells));   W % mean=0.
    allocate (P % mean(grid % n_cells));   P % mean=0.
  end if

  allocate (phi_face(grid % n_faces)); phi_face=0.

  allocate (phi_max(-grid % n_bnd_cells:grid % n_cells)); phi_max=0.
  allocate (phi_min(-grid % n_bnd_cells:grid % n_cells)); phi_min=0.

  allocate (G(6,grid % n_cells)); G=0

  allocate (Flux(grid % n_faces));     Flux=0.

  allocate (Utau(grid % n_materials));   Utau=0.0
  allocate (Vtau(grid % n_materials));   Vtau=0.0
  allocate (Wtau(grid % n_materials));   Wtau=0.0

  allocate (BadForG(grid % n_cells));  BadForG = .false.
  allocate (NumGood(grid % n_cells));  NumGood = 0          
  allocate (NumNeig(grid % n_cells));  NumNeig = 0         

  allocate (near(-grid % n_bnd_cells:grid % n_cells));  near  = 0.
  allocate (VISwall(-grid % n_bnd_cells:grid % n_cells)); VISwall =0.0

  ! For solution of temperature
  if(HOT==YES) then
    call Var_Mod_Allocate_Solution("T", t, grid)
    allocate (ut % n(-grid % n_bnd_cells:grid % n_cells));   ut % n=0.
    allocate (vt % n(-grid % n_bnd_cells:grid % n_cells));   vt % n=0.
    allocate (wt % n(-grid % n_bnd_cells:grid % n_cells));   wt % n=0.
    allocate (CONwall(-grid % n_bnd_cells:grid % n_cells)); CONwall =0.0
    if(BUOY == YES) then
      allocate (Pbuoy(-grid % n_bnd_cells:grid % n_cells));  Pbuoy=0.
    end if
  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(SIMULA==EBM.or.SIMULA==HJ) then
    if(SIMULA == HJ) then
      allocate (Eps_tot(-grid % n_bnd_cells:grid % n_cells)); Eps_tot=0.
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
  end if  ! SIMULA == EBM or HJ

  ! Variables for Rans models
  if(SIMULA==K_EPS.or.SIMULA == HYB_PITM) then

    call Var_Mod_Allocate_Solution("KIN", kin, grid)
    call Var_Mod_Allocate_Solution("EPS", eps, grid)

    allocate (Uf(-grid % n_bnd_cells:grid % n_cells));      Uf    =0.0
    allocate (Ufmean(-grid % n_bnd_cells:grid % n_cells));  Ufmean=0.0
    allocate (Pk(-grid % n_bnd_cells:grid % n_cells));      Pk    =0.0
    allocate (Ynd(-grid % n_bnd_cells:grid % n_cells));     Ynd   =0.0

    if(URANS == YES.or.SIMULA == HYB_PITM.or.SIMULA==HYB_ZETA) then
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

    if(URANS == YES.or.SIMULA==HYB_ZETA) then
      allocate (f22 % mean(grid % n_cells));   f22 % mean=0.
      allocate (v_2 % mean(grid % n_cells));   v_2 % mean=0.
    end if
  end if                    

  if(SIMULA == HYB_ZETA) then
    allocate (VISt_sgs(-grid % n_bnd_cells:grid % n_cells));  VISt_sgs=0.
    allocate (VISt_eff(-grid % n_bnd_cells:grid % n_cells));  VISt_eff=0.
  end if

  if(SIMULA == SPA_ALL.or.SIMULA == DES_SPA) then
    call Var_Mod_Allocate_Solution("VIS", vis, grid)
    if(SIMULA == DES_SPA) then
      allocate (Ksgs(-grid % n_bnd_cells:grid % n_cells));  Ksgs=0.
      allocate (VIS % mean(grid % n_cells));   VIS % mean=0.
    end if
  end if

  ! Variables defined in les_mod.h90:
  if(SIMULA == LES.or.SIMULA==HYB_ZETA) then
    if(MODE == WALE) then 
      allocate (WALEv(-grid % n_bnd_cells:grid % n_cells));  WALEv =0.
    end if
    if(MODE == DYN) then 
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

  if(SIMULA==LES.or.SIMULA==DNS.or.SIMULA==DES_SPA.or.&
     SIMULA==HYB_ZETA.or.SIMULA==HYB_PITM.or.URANS==YES) then
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
