!======================================================================!
  SUBROUTINE SavRes(namAut)
!----------------------------------------------------------------------!
! Writes restart files. NAME.restart                                                 !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, s, m
  CHARACTER           :: namOut*80, answer*80
  CHARACTER, OPTIONAL :: namAut*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: SavRes.f90,v 1.4 2017/09/06 13:13:27 mhadziabdic Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/SavRes.f90,v $  
!======================================================================!

  if(PRESENT(namAut)) then
!---- save the name
    answer = name
    name = namAut
  else
    if(this  < 2)                                                     &
      write(*,*) '# Output restart file name [skip cancels]:'
    call ReadC(7,inp,tn,ts,te)
!->>> write(*,*) inp(1:300)
    read(inp(ts(1):te(1)), '(A80)')  namOut
    answer=namOut
    call ToUppr(answer) 

    if(answer == 'SKIP') return 

!---- save the name
    answer = name
    name = namOut
  end if

!-----------------------------!
!     Create restart file     !
!-----------------------------!
  call NamFil(this, namOut, '.restart', len_trim('.restart') )
  open(9, FILE=namOut, FORM='UNFORMATTED')
  if(this  < 2) write(6, *) '# Now creating the file:', namOut

!---- version
  write(9) 0.0  ! version

!---- 60 INTEGER parameters -----------------------------------------
  write(9)        0,      NbC,       NC,       NS,     Ndtt,    Nstat
  write(9)       Cm,        0,        0,        0,        0,        0
  write(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  write(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,        0
  write(9)    BLEND,      HYB,   PER_BC,     MODE,     PIPE,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
!--------------------------------------------------------------------

!---- 60 REAL parameters --------------------------------------
  write(9)     0.0,    0.0,    0.0,     xp,     yp,     zp  
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0     
  write(9)   ReTau,   Tref,    Cs0,    0.0,    0.0,    0.0
  write(9)      dt,   Time,  Kflow,    0.0,    0.0,    0.0
  write(9)   U%URF,  P%URF,   URFC, SIMTol, U%Stol,    0.0 
  write(9) PP%Stol,  T%URF, T%STol,  Tflux,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
!--------------------------------------------------------------

  write(9) (U % n(c),   c=-NbC,NC)
  write(9) (V % n(c),   c=-NbC,NC)
  write(9) (W % n(c),   c=-NbC,NC)
  write(9) (U % o(c),   c=1,NC)
  write(9) (V % o(c),   c=1,NC)
  write(9) (W % o(c),   c=1,NC)

  write(9) (U % C(c),   c=1,NC)
  write(9) (V % C(c),   c=1,NC)
  write(9) (W % C(c),   c=1,NC)
  write(9) (U % Co(c),  c=1,NC)
  write(9) (V % Co(c),  c=1,NC)
  write(9) (W % Co(c),  c=1,NC)

  write(9) (U % Do(c),  c=1,NC)
  write(9) (V % Do(c),  c=1,NC)
  write(9) (W % Do(c),  c=1,NC)

  write(9) (U % X(c),   c=1,NC)
  write(9) (V % X(c),   c=1,NC)
  write(9) (W % X(c),   c=1,NC)
  write(9) (U % Xo(c),  c=1,NC)
  write(9) (V % Xo(c),  c=1,NC)
  write(9) (W % Xo(c),  c=1,NC)

  write(9) (P % n(c),   c=-NbC,NC)
  write(9) (PP % n(c),  c=-NbC,NC)

  write(9) (Px(c),   c=-NbC,NC)
  write(9) (Py(c),   c=-NbC,NC)
  write(9) (Pz(c),   c=-NbC,NC)

!---- Pressure drops in each material (domain)
  do m=1,Nmat
    write(9) PdropX(m), PdropY(m), PdropZ(m)
    write(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    write(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    write(9) AreaX(m),  AreaY(m),  AreaZ(m)
    write(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do


!---- Fluxes 
  write(9) (Flux(s), s=1,NS)

  if(HOT == YES) then
    write(9) (T % n(c),    c=-NbC,NC)
    write(9) (T % q(c),    c=-NbC,-1)
    write(9) (T % o(c),    c=1,NC)
    write(9) (T % C(c),    c=1,NC)
    write(9) (T % Co(c),   c=1,NC)
    write(9) (T % Do(c),   c=1,NC)
    write(9) (T % X(c),    c=1,NC)
    write(9) (T % Xo(c),   c=1,NC)
  end if

  if(SIMULA==K_EPS.or.SIMULA==ZETA.or.&
     SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA.or.SIMULA == HYB_PITM) then 
    write(9) (Kin % n(c),    c=-NbC,NC)
    write(9) (Kin % o(c),    c=1,NC)
    write(9) (Kin % C(c),    c=1,NC)
    write(9) (Kin % Co(c),   c=1,NC)
    write(9) (Kin % Do(c),   c=1,NC)
    write(9) (Kin % X(c),    c=1,NC)
    write(9) (Kin % Xo(c),   c=1,NC)

    write(9) (Eps % n(c),    c=-NbC,NC)
    write(9) (Eps % o(c),    c=1,NC)
    write(9) (Eps % C(c),    c=1,NC)
    write(9) (Eps % Co(c),   c=1,NC)
    write(9) (Eps % Do(c),   c=1,NC)
    write(9) (Eps % X(c),    c=1,NC)
    write(9) (Eps % Xo(c),   c=1,NC)

    write(9) (Pk(c),     c=-NbC,NC)
    write(9) (Uf(c),     c=-NbC,NC)
    write(9) (Ynd(c),    c=-NbC,NC) 
    write(9) (VISwall(c),c=-NbC,NC)
    write(9) (TauWall(c),c= 1,NC)
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA == HYB_ZETA) then
    write(9) (v_2%n(c),    c=-NbC,NC)
    write(9) (v_2%o(c),   c=1,NC)
    write(9) (v_2%C(c),   c=1,NC)
    write(9) (v_2%Co(c),  c=1,NC)
    write(9) (v_2%Do(c),  c=1,NC)
    write(9) (v_2%X(c),   c=1,NC)
    write(9) (v_2%Xo(c),  c=1,NC)

    write(9) (f22%n(c),    c=-NbC,NC)
    write(9) (f22%o(c),   c=1,NC)
    write(9) (f22%Do(c),  c=1,NC)
    write(9) (f22%X(c),   c=1,NC)
    write(9) (f22%Xo(c),  c=1,NC)
 
    write(9) (Tsc(c),    c=-NbC,NC)
    write(9) (Lsc(c),    c=-NbC,NC)
  end if 

  if(SIMULA==EBM.or.SIMULA==HJ) then
    write(9) (uu  % n(c),    c=-NbC,NC)
    write(9) (uu  % o(c),    c=1,NC)
    write(9) (uu  % C(c),    c=1,NC)
    write(9) (uu  % Co(c),   c=1,NC)
    write(9) (uu  % Do(c),   c=1,NC)
    write(9) (uu  % X(c),    c=1,NC)
    write(9) (uu  % Xo(c),   c=1,NC)

    write(9) (vv  % n(c),    c=-NbC,NC)
    write(9) (vv  % o(c),    c=1,NC)
    write(9) (vv  % C(c),    c=1,NC)
    write(9) (vv  % Co(c),   c=1,NC)
    write(9) (vv  % Do(c),   c=1,NC)
    write(9) (vv  % X(c),    c=1,NC)
    write(9) (vv  % Xo(c),   c=1,NC)

    write(9) (ww  % n(c),    c=-NbC,NC)
    write(9) (ww  % o(c),    c=1,NC)
    write(9) (ww  % C(c),    c=1,NC)
    write(9) (ww  % Co(c),   c=1,NC)
    write(9) (ww  % Do(c),   c=1,NC)
    write(9) (ww  % X(c),    c=1,NC)
    write(9) (ww  % Xo(c),   c=1,NC)

    write(9) (uv  % n(c),    c=-NbC,NC)
    write(9) (uv  % o(c),    c=1,NC)
    write(9) (uv  % C(c),    c=1,NC)
    write(9) (uv  % Co(c),   c=1,NC)
    write(9) (uv  % Do(c),   c=1,NC)
    write(9) (uv  % X(c),    c=1,NC)
    write(9) (uv  % Xo(c),   c=1,NC)

    write(9) (uw  % n(c),    c=-NbC,NC)
    write(9) (uw  % o(c),    c=1,NC)
    write(9) (uw  % C(c),    c=1,NC)
    write(9) (uw  % Co(c),   c=1,NC)
    write(9) (uw  % Do(c),   c=1,NC)
    write(9) (uw  % X(c),    c=1,NC)
    write(9) (uw  % Xo(c),   c=1,NC)

    write(9) (vw  % n(c),    c=-NbC,NC)
    write(9) (vw  % o(c),    c=1,NC)
    write(9) (vw  % C(c),    c=1,NC)
    write(9) (vw  % Co(c),   c=1,NC)
    write(9) (vw  % Do(c),   c=1,NC)
    write(9) (vw  % X(c),    c=1,NC)
    write(9) (vw  % Xo(c),   c=1,NC)

    write(9) (Eps % n(c),    c=-NbC,NC)
    write(9) (Eps % o(c),    c=1,NC)
    write(9) (Eps % C(c),    c=1,NC)
    write(9) (Eps % Co(c),   c=1,NC)
    write(9) (Eps % Do(c),   c=1,NC)
    write(9) (Eps % X(c),    c=1,NC)
    write(9) (Eps % Xo(c),   c=1,NC)

    write(9) (Pk(c),         c=-NbC,NC)
    write(9) (Kin % n(c),    c=-NbC,NC)
    write(9) (VISt(c),       c=-NbC,NC)

    if(SIMULA==EBM) then
      write(9) (f22 % n(c),    c=-NbC,NC)
      write(9) (f22 % o(c),    c=1,NC)
      write(9) (f22 % Do(c),   c=1,NC)
      write(9) (f22 % X(c),    c=1,NC)
      write(9) (f22 % Xo(c),   c=1,NC)
    end if

    if(URANS == YES) then
      write(9) (VAR10x(c),    c=-NbC,NC)
      write(9) (VAR10y(c),    c=-NbC,NC)
      write(9) (VAR10z(c),    c=-NbC,NC)
      write(9) (VAR11x(c),    c=-NbC,NC)
      write(9) (VAR11y(c),    c=-NbC,NC)
      write(9) (VAR11z(c),    c=-NbC,NC)
    end if  
  end if

 if(SIMULA == SPA_ALL.or.SIMULA==DES_SPA) then
    write(9) (VIS % n(c),    c=-NbC,NC)
    write(9) (VIS % o(c),    c=1,NC)
    write(9) (VIS % C(c),    c=1,NC)
    write(9) (VIS % Co(c),   c=1,NC)
    write(9) (VIS % Do(c),   c=1,NC)
    write(9) (VIS % X(c),    c=1,NC)
    write(9) (VIS % Xo(c),   c=1,NC)

    write(9) (Vort(c),       c=-NbC,NC)
  end if

  if(SIMULA/=DNS) write(9) (VISt(c), c=-NbC,NC)

  if(SIMULA==DNS.or.SIMULA==LES.or.SIMULA==URANS.or.&
     SIMULA==HYB_ZETA.or.SIMULA==HYB_PITM.or.SIMULA==DES_SPA) then
    write(9) (U % mean(c),  c=-NbC,NC)
    write(9) (V % mean(c),  c=-NbC,NC)
    write(9) (W % mean(c),  c=-NbC,NC)
    write(9) (uu % mean(c), c=-NbC,NC)
    write(9) (vv % mean(c), c=-NbC,NC)
    write(9) (ww % mean(c), c=-NbC,NC)
    write(9) (uv % mean(c), c=-NbC,NC)
    write(9) (uw % mean(c), c=-NbC,NC)
    write(9) (vw % mean(c), c=-NbC,NC)
    write(9) (P % mean(c),  c=1,NC)

    if(HOT == YES) then
      write(9) (T % mean(c), c=-NbC,NC)
      write(9) (TT % mean(c), c=-NbC,NC)
      write(9) (uT % mean(c), c=-NbC,NC)
      write(9) (vT % mean(c), c=-NbC,NC)
      write(9) (wT % mean(c), c=-NbC,NC)
    end if
  end if

  if(SIMULA == LES) then
    write(9) (VISt_mean(c), c=1,NC)
    write(9) (near(c),   c=-NbC,NC)
    if(MODE == DYN) write(9) (Cdyn_mean(c), c=1,NC)
  end if
 
  if(SIMULA == HYB_PITM.or.SIMULA==DES_SPA) then
    write(9) (VISt_mean(c), c=1,NC)
  end if

  if(SIMULA == HYB_ZETA) then
    write(9) (VISt_mean(c), c=1,NC)
    write(9) (VISt_sgs(c), c=1,NC)
    write(9) (VISt_eff(c), c=1,NC)
  end if

  close(9)

!---- restore the name
  name = answer 

  END SUBROUTINE SavRes 
