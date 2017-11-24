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
  use Tokenizer_Mod
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
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1), '(A80)')  name_in
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

  read(9) (U % a(c),    c = 1, grid % n_cells)
  read(9) (V % a(c),    c = 1, grid % n_cells)
  read(9) (W % a(c),    c = 1, grid % n_cells)
  read(9) (U % a_o(c),  c = 1, grid % n_cells)
  read(9) (V % a_o(c),  c = 1, grid % n_cells)
  read(9) (W % a_o(c),  c = 1, grid % n_cells)

  read(9) (U % d_o(c),  c = 1, grid % n_cells)
  read(9) (V % d_o(c),  c = 1, grid % n_cells)
  read(9) (W % d_o(c),  c = 1, grid % n_cells)

  read(9) (U % c(c),    c = 1, grid % n_cells)
  read(9) (V % c(c),    c = 1, grid % n_cells)
  read(9) (W % c(c),    c = 1, grid % n_cells)
  read(9) (U % c_o(c),  c = 1, grid % n_cells)
  read(9) (V % c_o(c),  c = 1, grid % n_cells)
  read(9) (W % c_o(c),  c = 1, grid % n_cells)

  read(9) (P % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (PP % n(c),  c = -grid % n_bnd_cells,grid % n_cells)

  read(9) (p % x(c),  c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (p % y(c),  c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (p % z(c),  c = -grid % n_bnd_cells,grid % n_cells)

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
    read(9) (T % a(c),    c = 1, grid % n_cells)
    read(9) (T % a_o(c),  c = 1, grid % n_cells)
    read(9) (T % d_o(c),  c = 1, grid % n_cells)
    read(9) (T % c(c),    c = 1, grid % n_cells)
    read(9) (T % c_o(c),  c = 1, grid % n_cells)
  end if

  if(SIMULA==K_EPS.or.SIMULA==ZETA.or.&
     SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA.or.SIMULA == HYB_PITM) then 
    read(9) (Kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Kin % o(c),    c = 1, grid % n_cells)
    read(9) (Kin % a(c),    c = 1, grid % n_cells)
    read(9) (Kin % a_o(c),  c = 1, grid % n_cells)
    read(9) (Kin % d_o(c),  c = 1, grid % n_cells)
    read(9) (Kin % c(c),    c = 1, grid % n_cells)
    read(9) (Kin % c_o(c),  c = 1, grid % n_cells)

    read(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Eps % o(c),    c = 1, grid % n_cells)
    read(9) (Eps % a(c),    c = 1, grid % n_cells)
    read(9) (Eps % a_o(c),  c = 1, grid % n_cells)
    read(9) (Eps % d_o(c),  c = 1, grid % n_cells)
    read(9) (Eps % c(c),    c = 1, grid % n_cells)
    read(9) (Eps % c_o(c),  c = 1, grid % n_cells)

    read(9) (Pk(c),       c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Uf(c),       c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Ynd(c),      c = -grid % n_bnd_cells,grid % n_cells) 
    read(9) (VISwall(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (TauWall(c),  c= 1,grid % n_cells)
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA == HYB_ZETA) then
    read(9) (v_2 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (v_2 % o(c),    c = 1, grid % n_cells)
    read(9) (v_2 % a(c),    c = 1, grid % n_cells)
    read(9) (v_2 % a_o(c),  c = 1, grid % n_cells)
    read(9) (v_2 % d_o(c),  c = 1, grid % n_cells)
    read(9) (v_2 % c(c),    c = 1, grid % n_cells)
    read(9) (v_2 % c_o(c),  c = 1, grid % n_cells)

    read(9) (f22 % n(c),     c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (f22 % o(c),     c = 1, grid % n_cells)
    read(9) (f22 % d_o(c),  c = 1, grid % n_cells)
    read(9) (f22 % c(c),    c = 1, grid % n_cells)
    read(9) (f22 % c_o(c),  c = 1, grid % n_cells)
 
    read(9) (Tsc(c), c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Lsc(c), c = -grid % n_bnd_cells,grid % n_cells)
  end if 

  if(SIMULA==EBM.or.SIMULA==HJ) then
    read(9) (uu % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uu % o(c),    c = 1, grid % n_cells)
    read(9) (uu % a(c),    c = 1, grid % n_cells)
    read(9) (uu % a_o(c),  c = 1, grid % n_cells)
    read(9) (uu % d_o(c),  c = 1, grid % n_cells)
    read(9) (uu % c(c),    c = 1, grid % n_cells)
    read(9) (uu % c_o(c),  c = 1, grid % n_cells)

    read(9) (vv % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vv % o(c),    c = 1, grid % n_cells)
    read(9) (vv % a(c),    c = 1, grid % n_cells)
    read(9) (vv % a_o(c),  c = 1, grid % n_cells)
    read(9) (vv % d_o(c),  c = 1, grid % n_cells)
    read(9) (vv % c(c),    c = 1, grid % n_cells)
    read(9) (vv % c_o(c),  c = 1, grid % n_cells)

    read(9) (ww % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (ww % o(c),    c = 1, grid % n_cells)
    read(9) (ww % a(c),    c = 1, grid % n_cells)
    read(9) (ww % a_o(c),  c = 1, grid % n_cells)
    read(9) (ww % d_o(c),  c = 1, grid % n_cells)
    read(9) (ww % c(c),    c = 1, grid % n_cells)
    read(9) (ww % c_o(c),  c = 1, grid % n_cells)

    read(9) (uv % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uv % o(c),    c = 1, grid % n_cells)
    read(9) (uv % a(c),    c = 1, grid % n_cells)
    read(9) (uv % a_o(c),  c = 1, grid % n_cells)
    read(9) (uv % d_o(c),  c = 1, grid % n_cells)
    read(9) (uv % c(c),    c = 1, grid % n_cells)
    read(9) (uv % c_o(c),  c = 1, grid % n_cells)

    read(9) (uw % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uw % o(c),    c = 1, grid % n_cells)
    read(9) (uw % a(c),    c = 1, grid % n_cells)
    read(9) (uw % a_o(c),  c = 1, grid % n_cells)
    read(9) (uw % d_o(c),  c = 1, grid % n_cells)
    read(9) (uw % c(c),    c = 1, grid % n_cells)
    read(9) (uw % c_o(c),  c = 1, grid % n_cells)

    read(9) (vw % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vw % o(c),    c = 1, grid % n_cells)
    read(9) (vw % a(c),    c = 1, grid % n_cells)
    read(9) (vw % a_o(c),  c = 1, grid % n_cells)
    read(9) (vw % d_o(c),  c = 1, grid % n_cells)
    read(9) (vw % c(c),    c = 1, grid % n_cells)
    read(9) (vw % c_o(c),  c = 1, grid % n_cells)

    read(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Eps % o(c),    c = 1, grid % n_cells)
    read(9) (Eps % a(c),    c = 1, grid % n_cells)
    read(9) (Eps % a_o(c),  c = 1, grid % n_cells)
    read(9) (Eps % d_o(c),  c = 1, grid % n_cells)
    read(9) (Eps % c(c),    c = 1, grid % n_cells)
    read(9) (Eps % c_o(c),   c = 1, grid % n_cells)

    read(9) (Pk(c),       c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Kin % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (VISt(c),     c = -grid % n_bnd_cells,grid % n_cells)

    if(SIMULA==EBM) then
      read(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (f22 % o(c),    c = 1, grid % n_cells)
      read(9) (f22 % d_o(c),  c = 1, grid % n_cells)
      read(9) (f22 % c(c),    c = 1, grid % n_cells)
      read(9) (f22 % c_o(c),  c = 1, grid % n_cells)
    end if

    if(URANS == YES) then
      read(9) (VAR10x(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR10y(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR10z(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11x(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11y(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (VAR11z(c),  c = -grid % n_bnd_cells,grid % n_cells)
    end if  
  end if

  if(SIMULA == SPA_ALL.or.SIMULA==DES_SPA) then
    read(9) (VIS % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (VIS % o(c),    c = 1, grid % n_cells)
    read(9) (VIS % a(c),    c = 1, grid % n_cells)
    read(9) (VIS % a_o(c),  c = 1, grid % n_cells)
    read(9) (VIS % d_o(c),  c = 1, grid % n_cells)
    read(9) (VIS % c(c),    c = 1, grid % n_cells)
    read(9) (VIS % c_o(c),  c = 1, grid % n_cells)

    read(9) (Vort(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(SIMULA/=DNS) read(9) (VISt(c), c = -grid % n_bnd_cells,grid % n_cells)

  if(SIMULA==DNS.or.SIMULA==LES.or.SIMULA==URANS.or.&
     SIMULA==HYB_ZETA.or.SIMULA==HYB_PITM.or.SIMULA==DES_SPA) then
    read(9) (U % mean(c),   c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (V % mean(c),   c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (W % mean(c),   c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uu % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vv % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (ww % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uv % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (uw % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vw % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (P % mean(c),   c = 1, grid % n_cells)

    if(HOT == YES) then
      read(9) (T % mean(c),   c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (TT % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (uT % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (vT % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (wT % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    end if
  end if

  if(SIMULA == LES) then
    read(9) (VISt_mean(c),  c = 1, grid % n_cells)
    read(9) (near(c),       c = -grid % n_bnd_cells,grid % n_cells)
    if(MODE == DYN) read(9) (Cdyn_mean(c), c = 1, grid % n_cells)
  end if
 
  if(SIMULA == HYB_PITM.or.SIMULA==DES_SPA) then
    read(9) (VISt_mean(c), c = 1, grid % n_cells)
  end if

  if(SIMULA == HYB_ZETA) then
    read(9) (VISt_mean(c), c = 1, grid % n_cells)
    read(9) (VISt_sgs(c),  c = 1, grid % n_cells)
    read(9) (VISt_eff(c),  c = 1, grid % n_cells)
  end if

  close(9)

  restart = .TRUE.

  ! Restore the name
  name = answer 

  end subroutine Load_Restart 
