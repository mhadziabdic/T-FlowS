!==============================================================================!
  subroutine Save_Restart(grid, name_aut)
!------------------------------------------------------------------------------!
!   Writes restart files. name.restart                                         !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod, only: this_proc
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Constants_Pro_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)     :: grid
  character, optional :: name_aut*(*)
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c, s, m
  character(len=80)   :: name_out, answer
!==============================================================================!

  if(PRESENT(name_aut)) then
    ! Save the name
    answer = problem_name
    problem_name = name_aut
  else
    if(this_proc  < 2)                                                     &
      print *, '# Output restart file name [skip cancels]:'
    call Tokenizer_Mod_Read_Line(CMN_FILE)
    read(line % tokens(1), '(A80)')  name_out
    answer=name_out
    call To_Upper_Case(answer) 

    if(answer == 'SKIP') return 

    ! save the name
    answer = problem_name
    problem_name = name_out
  end if

  !-------------------------!
  !   Create restart file   !
  !-------------------------!
  call Name_File(this_proc, name_out, '.restart')
  open(9, file=name_out, FORM='unformatted')
  if(this_proc  < 2) print *, '# Creating the file:', trim(name_out)

  ! Version
  write(9) 0.0  ! version

  ! 60 integer parameters 
  write(9)        0,      grid % n_bnd_cells,       grid % n_cells,       grid % n_faces,     Ndtt,    Nstat
  write(9)       Cm,        0,        0,        0,        0,        0
  write(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  write(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,        0
  write(9)    BLEND,        0,   PER_BC,     MODE,     PIPE,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0

  ! 60 real parameters
  write(9)     0.0,    0.0,    0.0,   bulk(1) % xp,  bulk(1) % yp,  bulk(1) % zp  
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0     
  write(9)   ReTau,   Tref,    Cs0,   Tinf,    0.0,    0.0
  write(9)      dt,   Time,  Kflow,    0.0,    0.0,    0.0
  write(9)   U%URF,  P%URF,   URFC, SIMTol, U%Stol,    0.0 
  write(9) PP%Stol,  T%URF, T%STol,  Tflux,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0

  write(9) (U % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (V % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (W % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (U % o(c),  c = 1, grid % n_cells)
  write(9) (V % o(c),  c = 1, grid % n_cells)
  write(9) (W % o(c),  c = 1, grid % n_cells)

  write(9) (U % a(c),    c = 1, grid % n_cells)
  write(9) (V % a(c),    c = 1, grid % n_cells)
  write(9) (W % a(c),    c = 1, grid % n_cells)
  write(9) (U % a_o(c),  c = 1, grid % n_cells)
  write(9) (V % a_o(c),  c = 1, grid % n_cells)
  write(9) (W % a_o(c),  c = 1, grid % n_cells)

  write(9) (U % d_o(c),  c = 1, grid % n_cells)
  write(9) (V % d_o(c),  c = 1, grid % n_cells)
  write(9) (W % d_o(c),  c = 1, grid % n_cells)

  write(9) (U % c(c),    c = 1, grid % n_cells)
  write(9) (V % c(c),    c = 1, grid % n_cells)
  write(9) (W % c(c),    c = 1, grid % n_cells)
  write(9) (U % c_o(c),  c = 1, grid % n_cells)
  write(9) (V % c_o(c),  c = 1, grid % n_cells)
  write(9) (W % c_o(c),  c = 1, grid % n_cells)

  write(9) (P % n(c),   c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (PP % n(c),  c = -grid % n_bnd_cells,grid % n_cells)

  write(9) (p % x(c),  c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (p % y(c),  c = -grid % n_bnd_cells,grid % n_cells)
  write(9) (p % z(c),  c = -grid % n_bnd_cells,grid % n_cells)

  ! Pressure drops in each material (domain)
  do m=1,grid % n_materials
    write(9) bulk(m) % p_drop_x,  bulk(m) % p_drop_y,  bulk(m) % p_drop_z
    write(9) bulk(m) % flux_x_o,  bulk(m) % flux_y_o,  bulk(m) % flux_z_o
    write(9) bulk(m) % flux_x,    bulk(m) % flux_y,    bulk(m) % flux_z
    write(9) bulk(m) % area_x,    bulk(m) % area_y,    bulk(m) % area_z
    write(9) bulk(m) % u,         bulk(m) % v,         bulk(m) % w
  end do

  ! Fluxes 
  write(9) (Flux(s), s=1,grid % n_faces)

  if(HOT == YES) then
    write(9) (T % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (T % q(c),    c = -grid % n_bnd_cells,-1)
    write(9) (T % o(c),    c = 1, grid % n_cells)
    write(9) (T % a(c),    c = 1, grid % n_cells)
    write(9) (T % a_o(c),  c = 1, grid % n_cells)
    write(9) (T % d_o(c),  c = 1, grid % n_cells)
    write(9) (T % c(c),    c = 1, grid % n_cells)
    write(9) (T % c_o(c),  c = 1, grid % n_cells)
  end if

  if(SIMULA==K_EPS.or.SIMULA==ZETA.or.&
     SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA.or.SIMULA == HYB_PITM) then 
    write(9) (Kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Kin % o(c),    c = 1, grid % n_cells)
    write(9) (Kin % a(c),    c = 1, grid % n_cells)
    write(9) (Kin % a_o(c),  c = 1, grid % n_cells)
    write(9) (Kin % d_o(c),  c = 1, grid % n_cells)
    write(9) (Kin % c(c),    c = 1, grid % n_cells)
    write(9) (Kin % c_o(c),  c = 1, grid % n_cells)

    write(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Eps % o(c),    c = 1, grid % n_cells)
    write(9) (Eps % a(c),    c = 1, grid % n_cells)
    write(9) (Eps % a_o(c),  c = 1, grid % n_cells)
    write(9) (Eps % d_o(c),  c = 1, grid % n_cells)
    write(9) (Eps % c(c),    c = 1, grid % n_cells)
    write(9) (Eps % c_o(c),  c = 1, grid % n_cells)

    write(9) (Pk(c),       c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Uf(c),       c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Ynd(c),      c = -grid % n_bnd_cells,grid % n_cells) 
    write(9) (VISwall(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (TauWall(c),  c= 1,grid % n_cells)
  end if

  if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA == HYB_ZETA) then
    write(9) (v_2 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (v_2 % o(c),    c = 1, grid % n_cells)
    write(9) (v_2 % a(c),    c = 1, grid % n_cells)
    write(9) (v_2 % a_o(c),  c = 1, grid % n_cells)
    write(9) (v_2 % d_o(c),  c = 1, grid % n_cells)
    write(9) (v_2 % c(c),    c = 1, grid % n_cells)
    write(9) (v_2 % c_o(c),  c = 1, grid % n_cells)

    write(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (f22 % o(c),    c = 1, grid % n_cells)
    write(9) (f22 % d_o(c),  c = 1, grid % n_cells)
    write(9) (f22 % c(c),    c = 1, grid % n_cells)
    write(9) (f22 % c_o(c),  c = 1, grid % n_cells)
 
    write(9) (Tsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Lsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if 

  if(SIMULA==EBM.or.SIMULA==HJ) then
    write(9) (uu % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uu % o(c),    c = 1, grid % n_cells)
    write(9) (uu % a(c),    c = 1, grid % n_cells)
    write(9) (uu % a_o(c),  c = 1, grid % n_cells)
    write(9) (uu % d_o(c),  c = 1, grid % n_cells)
    write(9) (uu % c(c),    c = 1, grid % n_cells)
    write(9) (uu % c_o(c),  c = 1, grid % n_cells)

    write(9) (vv % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vv % o(c),    c = 1, grid % n_cells)
    write(9) (vv % a(c),    c = 1, grid % n_cells)
    write(9) (vv % a_o(c),  c = 1, grid % n_cells)
    write(9) (vv % d_o(c),  c = 1, grid % n_cells)
    write(9) (vv % c(c),    c = 1, grid % n_cells)
    write(9) (vv % c_o(c),  c = 1, grid % n_cells)

    write(9) (ww % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (ww % o(c),    c = 1, grid % n_cells)
    write(9) (ww % a(c),    c = 1, grid % n_cells)
    write(9) (ww % a_o(c),  c = 1, grid % n_cells)
    write(9) (ww % d_o(c),  c = 1, grid % n_cells)
    write(9) (ww % c(c),    c = 1, grid % n_cells)
    write(9) (ww % c_o(c),  c = 1, grid % n_cells)

    write(9) (uv % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uv % o(c),    c = 1, grid % n_cells)
    write(9) (uv % a(c),    c = 1, grid % n_cells)
    write(9) (uv % a_o(c),  c = 1, grid % n_cells)
    write(9) (uv % d_o(c),  c = 1, grid % n_cells)
    write(9) (uv % c(c),    c = 1, grid % n_cells)
    write(9) (uv % c_o(c),  c = 1, grid % n_cells)

    write(9) (uw % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uw % o(c),    c = 1, grid % n_cells)
    write(9) (uw % a(c),    c = 1, grid % n_cells)
    write(9) (uw % a_o(c),  c = 1, grid % n_cells)
    write(9) (uw % d_o(c),  c = 1, grid % n_cells)
    write(9) (uw % c(c),    c = 1, grid % n_cells)
    write(9) (uw % c_o(c),  c = 1, grid % n_cells)

    write(9) (vw % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vw % o(c),    c = 1, grid % n_cells)
    write(9) (vw % a(c),    c = 1, grid % n_cells)
    write(9) (vw % a_o(c),  c = 1, grid % n_cells)
    write(9) (vw % d_o(c),  c = 1, grid % n_cells)
    write(9) (vw % c(c),    c = 1, grid % n_cells)
    write(9) (vw % c_o(c),  c = 1, grid % n_cells)

    write(9) (Eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Eps % o(c),    c = 1, grid % n_cells)
    write(9) (Eps % a(c),    c = 1, grid % n_cells)
    write(9) (Eps % a_o(c),  c = 1, grid % n_cells)
    write(9) (Eps % d_o(c),  c = 1, grid % n_cells)
    write(9) (Eps % c(c),    c = 1, grid % n_cells)
    write(9) (Eps % c_o(c),  c = 1, grid % n_cells)

    write(9) (Pk(c),         c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (VISt(c),       c = -grid % n_bnd_cells,grid % n_cells)

    if(SIMULA==EBM) then
      write(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (f22 % o(c),    c = 1, grid % n_cells)
      write(9) (f22 % d_o(c),  c = 1, grid % n_cells)
      write(9) (f22 % c(c),    c = 1, grid % n_cells)
      write(9) (f22 % c_o(c),  c = 1, grid % n_cells)
    end if

    if(URANS == YES) then
!     write(9) (VAR10x(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR10y(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR10z(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11x(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11y(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11z(c),  c = -grid % n_bnd_cells,grid % n_cells)
    end if  
  end if

  if(SIMULA == SPA_ALL.or.SIMULA==DES_SPA) then
    write(9) (VIS % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (VIS % o(c),    c = 1, grid % n_cells)
    write(9) (VIS % a(c),    c = 1, grid % n_cells)
    write(9) (VIS % a_o(c),  c = 1, grid % n_cells)
    write(9) (VIS % d_o(c),  c = 1, grid % n_cells)
    write(9) (VIS % c(c),    c = 1, grid % n_cells)
    write(9) (VIS % c_o(c),  c = 1, grid % n_cells)

    write(9) (Vort(c),       c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(SIMULA/=DNS) write(9) (VISt(c), c = -grid % n_bnd_cells,grid % n_cells)

  if(SIMULA==DNS.or.SIMULA==LES.or.SIMULA==URANS.or.&
     SIMULA==HYB_ZETA.or.SIMULA==HYB_PITM.or.SIMULA==DES_SPA) then
    write(9) (U % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (V % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (W % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uu % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vv % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (ww % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uv % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (uw % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vw % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (P % mean(c),  c = 1, grid % n_cells)

    if(HOT == YES) then
      write(9) (T % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (TT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (uT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (vT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (wT % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    end if
  end if

  if(SIMULA == LES) then
    write(9) (VISt_mean(c), c = 1, grid % n_cells)
    write(9) (near(c),   c = -grid % n_bnd_cells,grid % n_cells)
    if(MODE == DYN) write(9) (Cdyn_mean(c), c = 1, grid % n_cells)
  end if
 
  if(SIMULA == HYB_PITM.or.SIMULA==DES_SPA) then
    write(9) (VISt_mean(c), c = 1, grid % n_cells)
  end if

  if(SIMULA == HYB_ZETA) then
    write(9) (VISt_mean(c), c = 1, grid % n_cells)
    write(9) (VISt_sgs(c), c = 1, grid % n_cells)
    write(9) (VISt_eff(c), c = 1, grid % n_cells)
  end if

  close(9)

  ! Restore the name
  problem_name = answer 

  end subroutine
