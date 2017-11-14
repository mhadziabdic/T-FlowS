!==============================================================================!
  subroutine Load_Restart(grid, restart)
!------------------------------------------------------------------------------!
! Reads restart files name.restart                                             !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod, only: this_proc
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restart 
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s, m
  integer           :: i_1, i_2, i_3, i_4, i_5, i_6
  character(len=80) :: name_in, answer
  real              :: version
  real              :: r_1, r_2, r_3, r_4, r_5, r_6
!==============================================================================!

  if(this_proc  < 2) &              
    write(*,*) '# Input restart file name [skip cancels]:'
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)), '(A80)')  name_in
  answer=name_in
  call To_Upper_Case(answer) 

  if(answer == 'SKIP') then
    restart = .false.
    return 
  end if

  ! Save the name
  answer = name
  name = name_in

  !-----------------------!
  !   Read restart file   !
  !-----------------------!
  call Name_File(this_proc, name_in, '.restart', len_trim('.restart') )
  open(9, file=name_in, FORM='unformatted')
  write(6, *) '# Now reading the file:', name_in

  ! Version
  read(9) version ! version

  ! 60 integer parameters
  read(9)      i_1,      grid % n_bnd_cells,       grid % n_cells,       grid % n_faces,     Ndtt,    Nstat
  read(9)       Cm,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  read(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,      i_6
  read(9)    BLEND,      HYB,   PER_BC,     MODE,     PIPE,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6

  ! 60 real parameters
  read(9)     r_1,    r_2,    r_3,     xp,     yp,     zp  
  read(9)     r_1,    r_2,    r_3,    r_4,    r_4,    r_6             
  read(9)   ReTau,   Tref,    Cs0,   Tinf,    r_4,    r_6 
  read(9)      dt,   Time,  Kflow,    r_4,    r_5,    r_6 
  read(9)   U%URF,  P%URF,   URFC, SIMTol, U%Stol,    r_6 
  read(9) PP%Stol,  T%URF, T%STol,  Tflux,    r_5,    r_6   
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6 
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6   
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6   

  call Allocate_Variables(grid)

  read(9) (U % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (V % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (W % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (U % o(c),   c = 1, grid % n_cells)
  read(9) (V % o(c),   c = 1, grid % n_cells)
  read(9) (W % o(c),   c = 1, grid % n_cells)

  read(9) (U % C(c),   c = 1, grid % n_cells)
  read(9) (V % C(c),   c = 1, grid % n_cells)
  read(9) (W % C(c),   c = 1, grid % n_cells)
  read(9) (U % Co(c),  c = 1, grid % n_cells)
  read(9) (V % Co(c),  c = 1, grid % n_cells)
  read(9) (W % Co(c),  c = 1, grid % n_cells)

  read(9) (U % Do(c),  c = 1, grid % n_cells)
  read(9) (V % Do(c),  c = 1, grid % n_cells)
  read(9) (W % Do(c),  c = 1, grid % n_cells)

  read(9) (U % X(c),   c = 1, grid % n_cells)
  read(9) (V % X(c),   c = 1, grid % n_cells)
  read(9) (W % X(c),   c = 1, grid % n_cells)
  read(9) (U % Xo(c),  c = 1, grid % n_cells)
  read(9) (V % Xo(c),  c = 1, grid % n_cells)
  read(9) (W % Xo(c),  c = 1, grid % n_cells)

  read(9) (P % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (PP % n(c),  c = -grid % n_bnd_cells,grid % n_cells)

  read(9) (Px(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (Py(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (Pz(c),   c = -grid % n_bnd_cells,grid % n_cells)

  ! Pressure drops in each material (domain)
  do m = 1, grid % n_materials
    read(9) PdropX(m), PdropY(m), PdropZ(m)
    read(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    read(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    read(9) AreaX(m),  AreaY(m),  AreaZ(m)
    read(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do

  ! Fluxes 
  read(9) (Flux(s), s = 1, grid % n_faces)

  if(HOT == YES) then
    read(9) (T % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (T % q(c),    c = -grid % n_bnd_cells,-1)
    read(9) (T % o(c),    c = 1, grid % n_cells)
    read(9) (T % C(c),    c = 1, grid % n_cells)
    read(9) (T % Co(c),   c = 1, grid % n_cells)
    read(9) (T % Do(c),   c = 1, grid % n_cells)
    read(9) (T % X(c),    c = 1, grid % n_cells)
    read(9) (T % Xo(c),   c = 1, grid % n_cells)
  end if

  if(SIMULA==K_EPS.or.SIMULA==ZETA.or.&
     SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA.or.SIMULA == HYB_PITM) then 
    read(9) (Kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Kin % o(c),    c = 1, grid % n_cells)
    read(9) (Kin % C(c),    c = 1, grid % n_cells)
    read(9) (Kin % Co(c),   c = 1, grid % n_cells)
    read(9) (Kin % Do(c),   c = 1, grid % n_cells)
    read(9) (Kin % X(c),    c = 1, grid % n_cells)
    read(9) (Kin % Xo(c),   c = 1, grid % n_cells)

    read(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Eps % o(c),    c = 1, grid % n_cells)
    read(9) (Eps % C(c),    c = 1, grid % n_cells)
    read(9) (Eps % Co(c),   c = 1, grid % n_cells)
    read(9) (Eps % Do(c),   c = 1, grid % n_cells)
    read(9) (Eps % X(c),    c = 1, grid % n_cells)
    read(9) (Eps % Xo(c),   c = 1, grid % n_cells)

    read(9) (Pk(c),     c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Uf(c),     c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Ynd(c),    c = -grid % n_bnd_cells,grid % n_cells) 
    read(9) (VISwall(c),c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (TauWall(c),c= 1,grid % n_cells)
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA == HYB_ZETA) then
    read(9) (v_2%n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (v_2%o(c),   c = 1, grid % n_cells)
    read(9) (v_2%C(c),   c = 1, grid % n_cells)
    read(9) (v_2%Co(c),  c = 1, grid % n_cells)
    read(9) (v_2%Do(c),  c = 1, grid % n_cells)
    read(9) (v_2%X(c),   c = 1, grid % n_cells)
    read(9) (v_2%Xo(c),  c = 1, grid % n_cells)

    read(9) (f22%n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (f22%o(c),   c = 1, grid % n_cells)
    read(9) (f22%Do(c),  c = 1, grid % n_cells)
    read(9) (f22%X(c),   c = 1, grid % n_cells)
    read(9) (f22%Xo(c),  c = 1, grid % n_cells)
 
    read(9) (Tsc(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Lsc(c),    c = -grid % n_bnd_cells,grid % n_cells)
  end if 

  if(SIMULA==EBM.or.SIMULA==HJ) then
    read(9) (uu  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uu  % o(c),    c = 1, grid % n_cells)
    read(9) (uu  % C(c),    c = 1, grid % n_cells)
    read(9) (uu  % Co(c),   c = 1, grid % n_cells)
    read(9) (uu  % Do(c),   c = 1, grid % n_cells)
    read(9) (uu  % X(c),    c = 1, grid % n_cells)
    read(9) (uu  % Xo(c),   c = 1, grid % n_cells)

    read(9) (vv  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vv  % o(c),    c = 1, grid % n_cells)
    read(9) (vv  % C(c),    c = 1, grid % n_cells)
    read(9) (vv  % Co(c),   c = 1, grid % n_cells)
    read(9) (vv  % Do(c),   c = 1, grid % n_cells)
    read(9) (vv  % X(c),    c = 1, grid % n_cells)
    read(9) (vv  % Xo(c),   c = 1, grid % n_cells)

    read(9) (ww  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (ww  % o(c),    c = 1, grid % n_cells)
    read(9) (ww  % C(c),    c = 1, grid % n_cells)
    read(9) (ww  % Co(c),   c = 1, grid % n_cells)
    read(9) (ww  % Do(c),   c = 1, grid % n_cells)
    read(9) (ww  % X(c),    c = 1, grid % n_cells)
    read(9) (ww  % Xo(c),   c = 1, grid % n_cells)

    read(9) (uv  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uv  % o(c),    c = 1, grid % n_cells)
    read(9) (uv  % C(c),    c = 1, grid % n_cells)
    read(9) (uv  % Co(c),   c = 1, grid % n_cells)
    read(9) (uv  % Do(c),   c = 1, grid % n_cells)
    read(9) (uv  % X(c),    c = 1, grid % n_cells)
    read(9) (uv  % Xo(c),   c = 1, grid % n_cells)

    read(9) (uw  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uw  % o(c),    c = 1, grid % n_cells)
    read(9) (uw  % C(c),    c = 1, grid % n_cells)
    read(9) (uw  % Co(c),   c = 1, grid % n_cells)
    read(9) (uw  % Do(c),   c = 1, grid % n_cells)
    read(9) (uw  % X(c),    c = 1, grid % n_cells)
    read(9) (uw  % Xo(c),   c = 1, grid % n_cells)

    read(9) (vw  % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vw  % o(c),    c = 1, grid % n_cells)
    read(9) (vw  % C(c),    c = 1, grid % n_cells)
    read(9) (vw  % Co(c),   c = 1, grid % n_cells)
    read(9) (vw  % Do(c),   c = 1, grid % n_cells)
    read(9) (vw  % X(c),    c = 1, grid % n_cells)
    read(9) (vw  % Xo(c),   c = 1, grid % n_cells)

    read(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Eps % o(c),    c = 1, grid % n_cells)
    read(9) (Eps % C(c),    c = 1, grid % n_cells)
    read(9) (Eps % Co(c),   c = 1, grid % n_cells)
    read(9) (Eps % Do(c),   c = 1, grid % n_cells)
    read(9) (Eps % X(c),    c = 1, grid % n_cells)
    read(9) (Eps % Xo(c),   c = 1, grid % n_cells)

    read(9) (Pk(c),         c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (VISt(c),       c = -grid % n_bnd_cells,grid % n_cells)

    if(SIMULA==EBM) then
      read(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (f22 % o(c),    c = 1, grid % n_cells)
      read(9) (f22 % Do(c),   c = 1, grid % n_cells)
      read(9) (f22 % X(c),    c = 1, grid % n_cells)
      read(9) (f22 % Xo(c),   c = 1, grid % n_cells)
    end if

    if(URANS == YES) then
      read(9) (VAR10x(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR10y(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR10z(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11x(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11y(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11z(c),    c = -grid % n_bnd_cells,grid % n_cells)
    end if  
  end if

  if(SIMULA == SPA_ALL.or.SIMULA==DES_SPA) then
    read(9) (VIS % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (VIS % o(c),    c = 1, grid % n_cells)
    read(9) (VIS % C(c),    c = 1, grid % n_cells)
    read(9) (VIS % Co(c),   c = 1, grid % n_cells)
    read(9) (VIS % Do(c),   c = 1, grid % n_cells)
    read(9) (VIS % X(c),    c = 1, grid % n_cells)
    read(9) (VIS % Xo(c),   c = 1, grid % n_cells)

    read(9) (Vort(c),       c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(SIMULA/=DNS) read(9) (VISt(c), c = -grid % n_bnd_cells,grid % n_cells)

  if(SIMULA==DNS.or.SIMULA==LES.or.SIMULA==URANS.or.&
     SIMULA==HYB_ZETA.or.SIMULA==HYB_PITM.or.SIMULA==DES_SPA) then
    read(9) (U % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (V % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (W % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uu % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vv % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (ww % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uv % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uw % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vw % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (P % mean(c),  c = 1, grid % n_cells)

    if(HOT == YES) then
      read(9) (T % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (TT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (uT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (vT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (wT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    end if
  end if

  if(SIMULA == LES) then
    read(9) (VISt_mean(c), c = 1, grid % n_cells)
    read(9) (near(c),   c = -grid % n_bnd_cells,grid % n_cells)
    if(MODE == DYN) read(9) (Cdyn_mean(c), c = 1, grid % n_cells)
  end if
 
  if(SIMULA == HYB_PITM.or.SIMULA==DES_SPA) then
    read(9) (VISt_mean(c), c = 1, grid % n_cells)
  end if

  if(SIMULA == HYB_ZETA) then
    read(9) (VISt_mean(c), c = 1, grid % n_cells)
    read(9) (VISt_sgs(c), c = 1, grid % n_cells)
    read(9) (VISt_eff(c), c = 1, grid % n_cells)
  end if

  close(9)

  restart = .TRUE.

  ! Restore the name
  name = answer 

  end subroutine Load_Restart 
