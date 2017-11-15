!==============================================================================!
  subroutine Load_Restart_Ini(grid)
!------------------------------------------------------------------------------!
! Reads: name.restart                                                          !
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
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s, m
  integer           :: i_1, i_2, i_3, i_4, i_5, i_6
  character(len=80) :: name_in, answer
  real              :: version
  real              :: r_1, r_2, r_3, r_4, r_5, r_6
!==============================================================================!

  if(this_proc  < 2) &              
    write(*,*) '# Input intial restart file name [write skip to continue]:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(token % string(token % s(1):token % e(1)), '(A80)')  name_in
  answer=name_in
  call To_Upper_Case(answer) 

  if(answer == 'SKIP') then
    return 
  end if

  ! Initiated field from previous computation 
  if(this_proc  < 2) then
    write(*,*) '# Initialization of fields from previous computation: '
    write(*,*) '# DNS      -> Direct Numerical Simulation'
    write(*,*) '# LES      -> Large Eddy Simulation'
    write(*,*) '# K_EPS    -> High Reynolds k-eps model.'
    write(*,*) '# K_EPS_VV -> Durbin`s model.'
    write(*,*) '# SPA_ALL  -> Spalart-Allmaras model.'
    write(*,*) '# ZETA     -> k-eps-zeta-f model.'
    write(*,*) '# HJ       -> HJ model.'
    write(*,*) '# EBM      -> EBM model.'
  endif
  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(token % string(token % s(1):token % e(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'DNS') then
    RES_INI = DNS
  else if(answer == 'LES') then
    RES_INI = LES
  else if(answer == 'K_EPS') then
    RES_INI = K_EPS
  else if(answer == 'K_EPS_VV') then
    RES_INI = K_EPS_VV
  else if(answer == 'SPA_ALL') then
    RES_INI = SPA_ALL
  else if(answer == 'DES_SPA') then
    RES_INI = DES_SPA
  else if(answer == 'ZETA') then
    RES_INI = ZETA
  else if(answer == 'HYB_PITM') then
    RES_INI = HYB_PITM
  else if(answer == 'HYB_ZETA') then
    RES_INI = HYB_ZETA
  else if(answer == 'EBM') then
    RES_INI = EBM
  else if(answer == 'HJ') then
    RES_INI = HJ
  else if(answer == 'SKIP') then
    RES_INI = 1000 
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

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
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6

  ! 60 real parameters 
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6 
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6   
  read(9)     r_1,    r_2,    r_3,    r_4,    r_5,    r_6   

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
  do m=1,grid % n_materials
    read(9) PdropX(m), PdropY(m), PdropZ(m)
    read(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    read(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    read(9) AreaX(m),  AreaY(m),  AreaZ(m)
    read(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do

  ! Fluxes 
  read(9) (Flux(s), s=1,grid % n_faces)

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
  
  if(RES_INI==K_EPS.or.RES_INI==ZETA.or.&
     RES_INI==K_EPS_VV.or.RES_INI==HYB_ZETA.or.RES_INI == HYB_PITM) then 
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
    read(9) (TauWall(c),c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(RES_INI==K_EPS_VV.or.RES_INI==ZETA.or.RES_INI == HYB_ZETA) then
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

  if(RES_INI==EBM.or.RES_INI==HJ) then
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

    if(RES_INI==EBM) then
      read(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (f22 % o(c),    c = 1, grid % n_cells)
      read(9) (f22 % Do(c),   c = 1, grid % n_cells)
      read(9) (f22 % X(c),    c = 1, grid % n_cells)
      read(9) (f22 % Xo(c),   c = 1, grid % n_cells)
    end if
  end if

  if(RES_INI == SPA_ALL.or.RES_INI==DES_SPA) then
    read(9) (VIS % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (VIS % o(c),    c = 1, grid % n_cells)
    read(9) (VIS % C(c),    c = 1, grid % n_cells)
    read(9) (VIS % Co(c),   c = 1, grid % n_cells)
    read(9) (VIS % Do(c),   c = 1, grid % n_cells)
    read(9) (VIS % X(c),    c = 1, grid % n_cells)
    read(9) (VIS % Xo(c),   c = 1, grid % n_cells)

    read(9) (Vort(c),       c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(RES_INI/=DNS) read(9) (VISt(c), c = -grid % n_bnd_cells,grid % n_cells)

  close(9)

  ! Restore the name
  name = answer 

  write(*,*) 'Leaving Load_Restart_Ini.f90'

  end subroutine Load_Restart_Ini
