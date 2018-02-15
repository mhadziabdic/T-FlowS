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
  use Control_Mod
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
  character(len=80) :: old_heat_transfer
  character(len=80) :: old_turbulence_model
!==============================================================================!

  call Control_Mod_Load_Restart_Name(name_in)

  answer=name_in
  call To_Upper_Case(answer)
  if(answer == 'SKIP') return

  call Control_Mod_Heat_Transfer(old_heat_transfer)
  call Control_Mod_Turbulence_Model(old_turbulence_model)

  ! Save the name
  answer = problem_name
  problem_name = name_in

  !-----------------------!
  !   Read restart file   !
  !-----------------------!
  call Name_File(this_proc, name_in, '.restart')
  open(9, file=name_in, form='unformatted')
  print *, '# Reading the file: ', name_in

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

  read(9) (U % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (V % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (W % n(c),  c = -grid % n_bnd_cells,grid % n_cells)
  read(9) (U % o(c),  c = 1, grid % n_cells)
  read(9) (V % o(c),  c = 1, grid % n_cells)
  read(9) (W % o(c),  c = 1, grid % n_cells)

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
  do m=1,grid % n_materials
    read(9) bulk(m) % p_drop_x,  bulk(m) % p_drop_y,  bulk(m) % p_drop_z
    read(9) bulk(m) % flux_x_o,  bulk(m) % flux_y_o,  bulk(m) % flux_z_o
    read(9) bulk(m) % flux_x,    bulk(m) % flux_y,    bulk(m) % flux_z
    read(9) bulk(m) % area_x,    bulk(m) % area_y,    bulk(m) % area_z
    read(9) bulk(m) % u,         bulk(m) % v,         bulk(m) % w
  end do

  ! Fluxes 
  read(9) (Flux(s), s=1,grid % n_faces)

  if(old_heat_transfer == 'YES') then
    read(9) (T % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (T % q(c),    c = -grid % n_bnd_cells,-1)
    read(9) (T % o(c),    c = 1, grid % n_cells)
    read(9) (T % a(c),    c = 1, grid % n_cells)
    read(9) (T % a_o(c),  c = 1, grid % n_cells)
    read(9) (T % d_o(c),  c = 1, grid % n_cells)
    read(9) (T % c(c),    c = 1, grid % n_cells)
    read(9) (T % c_o(c),  c = 1, grid % n_cells)
  end if
  
  if(old_turbulence_model == 'K_EPS'    .or.  &
     old_turbulence_model == 'ZETA'     .or.  &
     old_turbulence_model == 'K_EPS_VV' .or.  &
     old_turbulence_model == 'HYB_ZETA' .or.  &
     old_turbulence_model == 'HYB_PITM') then 
    read(9) (kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (kin % o(c),    c = 1, grid % n_cells)
    read(9) (kin % a(c),    c = 1, grid % n_cells)
    read(9) (kin % a_o(c),  c = 1, grid % n_cells)
    read(9) (kin % d_o(c),  c = 1, grid % n_cells)
    read(9) (kin % c(c),    c = 1, grid % n_cells)
    read(9) (kin % c_o(c),  c = 1, grid % n_cells)

    read(9) (eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (eps % o(c),    c = 1, grid % n_cells)
    read(9) (eps % a(c),    c = 1, grid % n_cells)
    read(9) (eps % a_o(c),  c = 1, grid % n_cells)
    read(9) (eps % d_o(c),  c = 1, grid % n_cells)
    read(9) (eps % c(c),    c = 1, grid % n_cells)
    read(9) (eps % c_o(c),  c = 1, grid % n_cells)

    read(9) (p_kin(c),       c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Uf(c),       c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Ynd(c),      c = -grid % n_bnd_cells,grid % n_cells) 
    read(9) (viswall(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (TauWall(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(old_turbulence_model == 'K_EPS_VV' .or.  &
     old_turbulence_model == 'ZETA'     .or.  &
     old_turbulence_model == 'HYB_ZETA') then
    read(9) (v_2 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (v_2 % o(c),    c = 1, grid % n_cells)
    read(9) (v_2 % a(c),    c = 1, grid % n_cells)
    read(9) (v_2 % a_o(c),  c = 1, grid % n_cells)
    read(9) (v_2 % d_o(c),  c = 1, grid % n_cells)
    read(9) (v_2 % c(c),    c = 1, grid % n_cells)
    read(9) (v_2 % c_o(c),  c = 1, grid % n_cells)

    read(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (f22 % o(c),    c = 1, grid % n_cells)
    read(9) (f22 % d_o(c),  c = 1, grid % n_cells)
    read(9) (f22 % c(c),    c = 1, grid % n_cells)
    read(9) (f22 % c_o(c),  c = 1, grid % n_cells)
 
    read(9) (Tsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (Lsc(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if 

  if(old_turbulence_model == 'EBM' .or.  &
     old_turbulence_model == 'HJ') then
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

    read(9) (eps % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (eps % o(c),    c = 1, grid % n_cells)
    read(9) (eps % a(c),    c = 1, grid % n_cells)
    read(9) (eps % a_o(c),  c = 1, grid % n_cells)
    read(9) (eps % d_o(c),  c = 1, grid % n_cells)
    read(9) (eps % c(c),    c = 1, grid % n_cells)
    read(9) (eps % c_o(c),  c = 1, grid % n_cells)

    read(9) (p_kin(c),      c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (kin % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vis_t(c),       c = -grid % n_bnd_cells,grid % n_cells)

    if(old_turbulence_model == 'EBM') then
      read(9) (f22 % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
      read(9) (f22 % o(c),    c = 1, grid % n_cells)
      read(9) (f22 % d_o(c),  c = 1, grid % n_cells)
      read(9) (f22 % c(c),    c = 1, grid % n_cells)
      read(9) (f22 % c_o(c),  c = 1, grid % n_cells)
    end if
  end if

  if(old_turbulence_model == 'SPA_ALL' .or.  &
     old_turbulence_model == 'DES_SPA') then
    read(9) (vis % n(c),    c = -grid % n_bnd_cells,grid % n_cells)
    read(9) (vis % o(c),    c = 1, grid % n_cells)
    read(9) (vis % a(c),    c = 1, grid % n_cells)
    read(9) (vis % a_o(c),  c = 1, grid % n_cells)
    read(9) (vis % d_o(c),  c = 1, grid % n_cells)
    read(9) (vis % c(c),    c = 1, grid % n_cells)
    read(9) (vis % c_o(c),  c = 1, grid % n_cells)

    read(9) (vort(c),  c = -grid % n_bnd_cells,grid % n_cells)
  end if

  if(old_turbulence_model /= 'DNS')  &
    read(9) (vis_t(c), c = -grid % n_bnd_cells,grid % n_cells)

  close(9)

  ! Restore the name
  problem_name = answer 

  print *, 'Leaving Load_Restart_Ini.f90'

  end subroutine
