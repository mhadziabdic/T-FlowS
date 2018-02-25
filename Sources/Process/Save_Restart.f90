!==============================================================================!
  subroutine Save_Restart(grid, time_step, name_aut)
!------------------------------------------------------------------------------!
!   Writes restart files. name.restart                                         !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use les_mod
  use Comm_Mod, only: this_proc
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)     :: grid
  integer             :: time_step
  character, optional :: name_aut*(*)
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s, m
  character(len=80) :: name_out, answer
!==============================================================================!

  if(present(name_aut)) then
    ! Save the name
    answer = problem_name
    problem_name = name_aut
  else
    call Control_Mod_Load_Restart_Name(name_out)

    answer=name_out
    call To_Upper_Case(answer)
    if(answer == 'SKIP') then
      return 
    end if

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
  if(this_proc  < 2) print *, '# Creating the file: ', trim(name_out)

  ! Version
  write(9) 0.0  ! version

  ! 60 integer parameters 
  write(9) time_step,     grid % n_bnd_cells,       grid % n_cells,       grid % n_faces,     0,    0
  write(9)       Cm,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)   POSPRO,  CHANNEL,     TEST,    OTHER,        0,        0
  write(9)        0,        0,   PER_BC,        0,     PIPE,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0

  ! 60 real parameters
  write(9)     0.0,    0.0,    0.0,   bulk(1) % xp,  bulk(1) % yp,  bulk(1) % zp  
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0     
  write(9)   ReTau,   Tref,    Cs0,   Tinf,    0.0,    0.0
  write(9)     0.0,    0.0,  Kflow,    0.0,    0.0,    0.0
  write(9)     0.0,    0.0,    0.0,    0.0,    0.0,    0.0 
  write(9)     0.0,    0.0,    0.0,  Tflux,    0.0,    0.0
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
  write(9) (flux(s), s=1,grid % n_faces)

  if(heat_transfer == YES) then
    write(9) (T % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (T % q(c),    c = -grid % n_bnd_cells,-1)
    write(9) (T % o(c),    c = 1, grid % n_cells)
    write(9) (T % a(c),    c = 1, grid % n_cells)
    write(9) (T % a_o(c),  c = 1, grid % n_cells)
    write(9) (T % d_o(c),  c = 1, grid % n_cells)
    write(9) (T % c(c),    c = 1, grid % n_cells)
    write(9) (T % c_o(c),  c = 1, grid % n_cells)
  end if

  if(turbulence_model == K_EPS    .or.  &
     turbulence_model == K_EPS_ZETA_F     .or.  &
     turbulence_model == K_EPS_V2 .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_PITM) then 
    write(9) (kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (kin % o(c),    c = 1, grid % n_cells)
    write(9) (kin % a(c),    c = 1, grid % n_cells)
    write(9) (kin % a_o(c),  c = 1, grid % n_cells)
    write(9) (kin % d_o(c),  c = 1, grid % n_cells)
    write(9) (kin % c(c),    c = 1, grid % n_cells)
    write(9) (kin % c_o(c),  c = 1, grid % n_cells)

    write(9) (eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (eps % o(c),    c = 1, grid % n_cells)
    write(9) (eps % a(c),    c = 1, grid % n_cells)
    write(9) (eps % a_o(c),  c = 1, grid % n_cells)
    write(9) (eps % d_o(c),  c = 1, grid % n_cells)
    write(9) (eps % c(c),    c = 1, grid % n_cells)
    write(9) (eps % c_o(c),  c = 1, grid % n_cells)

    write(9) (p_kin(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Uf(c),       c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Ynd(c),      c = -grid % n_bnd_cells,grid % n_cells) 
    write(9) (vis_wall(c), c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (tau_wall(c), c = 1, grid % n_cells)
  end if

  if(turbulence_model == K_EPS_V2 .or.  &
     turbulence_model == K_EPS_ZETA_F     .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    write(9) (v2 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (v2 % o(c),    c = 1, grid % n_cells)
    write(9) (v2 % a(c),    c = 1, grid % n_cells)
    write(9) (v2 % a_o(c),  c = 1, grid % n_cells)
    write(9) (v2 % d_o(c),  c = 1, grid % n_cells)
    write(9) (v2 % c(c),    c = 1, grid % n_cells)
    write(9) (v2 % c_o(c),  c = 1, grid % n_cells)

    write(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (f22 % o(c),    c = 1, grid % n_cells)
    write(9) (f22 % d_o(c),  c = 1, grid % n_cells)
    write(9) (f22 % c(c),    c = 1, grid % n_cells)
    write(9) (f22 % c_o(c),  c = 1, grid % n_cells)
 
    write(9) (Tsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (Lsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if 

  if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC) then
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

    write(9) (eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (eps % o(c),    c = 1, grid % n_cells)
    write(9) (eps % a(c),    c = 1, grid % n_cells)
    write(9) (eps % a_o(c),  c = 1, grid % n_cells)
    write(9) (eps % d_o(c),  c = 1, grid % n_cells)
    write(9) (eps % c(c),    c = 1, grid % n_cells)
    write(9) (eps % c_o(c),  c = 1, grid % n_cells)

    write(9) (p_kin(c),      c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vis_t(c),      c = -grid % n_bnd_cells,grid % n_cells)

    if(turbulence_model == REYNOLDS_STRESS_MODEL) then
      write(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (f22 % o(c),    c = 1, grid % n_cells)
      write(9) (f22 % d_o(c),  c = 1, grid % n_cells)
      write(9) (f22 % c(c),    c = 1, grid % n_cells)
      write(9) (f22 % c_o(c),  c = 1, grid % n_cells)
    end if

    if(turbulence_model_variant == URANS) then
!     write(9) (VAR10x(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR10y(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR10z(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11x(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11y(c),  c = -grid % n_bnd_cells,grid % n_cells)
!     write(9) (VAR11z(c),  c = -grid % n_bnd_cells,grid % n_cells)
    end if  
  end if

  if(turbulence_model == SPALART_ALLMARAS .or.  &
     turbulence_model == DES_SPALART) then
    write(9) (vis % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    write(9) (vis % o(c),    c = 1, grid % n_cells)
    write(9) (vis % a(c),    c = 1, grid % n_cells)
    write(9) (vis % a_o(c),  c = 1, grid % n_cells)
    write(9) (vis % d_o(c),  c = 1, grid % n_cells)
    write(9) (vis % c(c),    c = 1, grid % n_cells)
    write(9) (vis % c_o(c),  c = 1, grid % n_cells)

    write(9) (Vort(c),       c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(turbulence_model /= DNS)  &
    write(9) (vis_t(c), c = -grid % n_bnd_cells,grid % n_cells)

  if(turbulence_model         == DNS                 .or.  &
     turbulence_model         == LES                 .or.  &
     turbulence_model         == HYBRID_K_EPS_ZETA_F .or.  &
     turbulence_model         == HYBRID_PITM         .or.  &
     turbulence_model         == DES_SPALART         .or.  &
     turbulence_model_variant == URANS) then
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

    if(heat_transfer == YES) then
      write(9) (t % mean(c),  c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (tt % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (ut % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (vt % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
      write(9) (wt % mean(c), c = -grid % n_bnd_cells,grid % n_cells)
    end if
  end if

  if(turbulence_model == LES) then
    write(9) (vis_t_mean(c),        c = 1, grid % n_cells)
    write(9) (nearest_wall_cell(c), c = 1, grid % n_cells)
    if(turbulence_model_variant == DYNAMIC)  &
      write(9) (c_dyn_mean(c), c = 1, grid % n_cells)
  end if
 
  if(turbulence_model == HYBRID_PITM .or.  &
     turbulence_model == DES_SPALART) then
    write(9) (vis_t_mean(c), c = 1, grid % n_cells)
  end if

  if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    write(9) (vis_t_mean(c), c = 1, grid % n_cells)
    write(9) (vis_t_sgs(c),  c = 1, grid % n_cells)
    write(9) (vis_t_eff(c),  c = 1, grid % n_cells)
  end if

  close(9)

  ! Restore the name
  problem_name = answer 

  end subroutine
