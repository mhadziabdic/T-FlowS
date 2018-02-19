!==============================================================================!
  subroutine Cgns_Mod_Add_Fields_To_Grid_Seq(grid, name_save)
!------------------------------------------------------------------------------!
!   Adds fields to existing grid cgns file [sequential vesion]                 !
!------------------------------------------------------------------------------!
  use allp_mod
  use all_mod
  use Flow_Mod
  use rans_mod
  use par_mod, only: this_proc, n_proc
  use Tokenizer_Mod
  use Grid_Mod
  use Cgns_Mod
  use Work_Mod, only: uu_mean => r_cell_01,  &
                      vv_mean => r_cell_02,  &
                      ww_mean => r_cell_03,  &
                      uv_mean => r_cell_04,  &
                      uw_mean => r_cell_05,  &
                      vw_mean => r_cell_06,  &
                      tt_mean => r_cell_07,  &
                      ut_mean => r_cell_08,  &
                      vt_mean => r_cell_09,  &
                      wt_mean => r_cell_10
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: store_name
  integer           :: base
  integer           :: block
  integer           :: solution
  integer           :: field
  integer           :: c
!==============================================================================!

  print *, "# subroutine Cgns_Mod_Save_Grid_Seq"

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  !--------------------------!
  !   Open file for modify   !
  !--------------------------!
  call Name_File(0, file_name, '.cgns')

  file_mode = CG_MODE_MODIFY
  call Cgns_Mod_Open_File_Seq(file_mode)

  call Initialize_Counters

  ! Count number of 3d cell type elements
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 8) cnt_hex = cnt_hex + 1
    if(grid % cells_n_nodes(c) == 6) cnt_wed = cnt_wed + 1
    if(grid % cells_n_nodes(c) == 5) cnt_pyr = cnt_pyr + 1
    if(grid % cells_n_nodes(c) == 4) cnt_tet = cnt_tet + 1
  end do

  !-----------------!
  !                 !
  !   Bases block   !
  !                 !
  !-----------------!
  n_bases = 1
  allocate(cgns_base(n_bases))

  base = 1
  cgns_base(base) % name = "Base 1"
  cgns_base(base) % cell_dim = 3
  cgns_base(base) % phys_dim = 3

  !-----------------!
  !                 !
  !   Zones block   !
  !                 !
  !-----------------!

  cgns_base(base) % n_blocks = 1
  allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

  block = 1
  cgns_base(base) % block(block) % name = "Zone 1"
  cgns_base(base) % block(block) % mesh_info(1) = grid % n_nodes
  cgns_base(base) % block(block) % mesh_info(2) = grid % n_cells
  cgns_base(base) % block(block) % mesh_info(3) = 0

  !--------------------!
  !                    !
  !   Solution block   !
  !                    !
  !--------------------!

  cgns_base(base) % block(block) % n_solutions = 1

  allocate(cgns_base(base) % block(block) % solution( &
    cgns_base(base) % block(block) % n_solutions))
  solution = 1

  cgns_base(base) % block(block) % solution(solution) % name = 'FlowSolution'
  cgns_base(base) % block(block) % solution(solution) % sol_type = CellCenter

  call Cgns_Mod_Write_Solution_Info_Seq(base, block, solution)

  !-----------------!
  !                 !
  !   Field block   !
  !                 !
  !-----------------!

  !-------------------------------------------!
  !   Copy code below from Save_Vtu_Results   !
  !-------------------------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    U % n(1:grid % n_cells),'VelocityX')
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    V % n(1:grid % n_cells),'VelocityY')
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    W % n(1:grid % n_cells),'VelocityZ')
  !--------------!
  !   Pressure   !
  !--------------!
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    P % n(1:grid % n_cells),'Pressure')
  !-----------------!
  !   Temperature   !
  !-----------------!
  if(heat_transfer == YES) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      T % n(1:grid % n_cells),'Temperature')
  end if
  !--------------------------!
  !   Turbulent quantities   !
  !--------------------------!

  ! Kin and Eps
  if(turbulence_model == K_EPS                 .or.  &
     turbulence_model == K_EPS_V2              .or.  &
     turbulence_model == K_EPS_ZETA_F          .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F   .or.  &
     turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC  ) then

    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      Kin % n(1:grid % n_cells),'TurbulentEnergyKinetic')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      Eps % n(1:grid % n_cells),'TurbulentDissipation')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      p_kin(1:grid % n_cells),'TurbulentEnergyKineticProduction')
  end if

  ! v2 and f22
  if(turbulence_model == K_EPS_V2       .or.  &
     turbulence_model == K_EPS_ZETA_F   .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      v2 % n(1:grid % n_cells), v2 % name)
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      f22 % n(1:grid % n_cells), f22 % name)
  end if

  ! Vis and Turbulent Vicosity_t
  if(turbulence_model == DES_SPALART .or.  &
     turbulence_model == SPALART_ALLMARAS) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vis % n(1:grid % n_cells),'ViscosityKinematic')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vort(1:grid % n_cells),'VorticityMagnitude')
  end if
  if(turbulence_model == K_EPS                 .or.  &
     turbulence_model == K_EPS_V2              .or.  &
     turbulence_model == K_EPS_ZETA_F          .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F   .or.  &
     turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC     .or.  &
     turbulence_model == LES                   .or.  &
     turbulence_model == DES_SPALART           .or.  &
     turbulence_model == SPALART_ALLMARAS) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vis_t(1:grid % n_cells),'ViscosityEddyKinematic')
  end if

  ! Reynolds stress models
  if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uu % n(1:grid % n_cells),'ReynoldsStressXX')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vv % n(1:grid % n_cells),'ReynoldsStressYY')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      ww % n(1:grid % n_cells),'ReynoldsStressZZ')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uv % n(1:grid % n_cells),'ReynoldsStressXY')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uw % n(1:grid % n_cells),'ReynoldsStressXZ')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vw % n(1:grid % n_cells),'ReynoldsStressYZ')
  end if

  ! Statistics for large-scale simulations of turbulence
  if(turbulence_model == LES .or.  &
     turbulence_model == DES_SPALART) then
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      U % n(1:grid % n_cells),'Velocity_MeanX')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      U % n(1:grid % n_cells),'Velocity_MeanY')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      U % n(1:grid % n_cells),'Velocity_MeanZ')
    uu_mean = uu % mean(c) - u % mean(c) * u % mean(c)
    vv_mean = vv % mean(c) - v % mean(c) * v % mean(c)
    ww_mean = ww % mean(c) - w % mean(c) * w % mean(c)
    uv_mean = uv % mean(c) - u % mean(c) * v % mean(c)
    uw_mean = uw % mean(c) - u % mean(c) * w % mean(c)
    vw_mean = vw % mean(c) - v % mean(c) * w % mean(c)
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uu_mean(1:grid % n_cells),'ReynoldsStress_MeanXX')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vv_mean(1:grid % n_cells),'ReynoldsStress_MeanYY')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      ww_mean(1:grid % n_cells),'ReynoldsStress_MeanZZ')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uv_mean(1:grid % n_cells),'ReynoldsStress_MeanXY')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      uw_mean(1:grid % n_cells),'ReynoldsStress_MeanXZ')
    call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
      vw_mean(1:grid % n_cells),'ReynoldsStress_MeanYZ')
    if(heat_transfer == YES) then
      call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
        t % mean(1:grid % n_cells),'Temperature_Mean')
      tt_mean = tt % mean(c) - t % mean(c) * t % mean(c)
      ut_mean = ut % mean(c) - u % mean(c) * t % mean(c)
      vt_mean = vt % mean(c) - v % mean(c) * t % mean(c)
      wt_mean = wt % mean(c) - w % mean(c) * t % mean(c)
      call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
        uu_mean(1:grid % n_cells),'TT_Mean')
      call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
        ut_mean(1:grid % n_cells),'UT_Mean')
      call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
        vt_mean(1:grid % n_cells),'VT_Mean')
      call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
        wt_mean(1:grid % n_cells),'WT_Mean')
    end if
  end if

  ! Wall distance and delta
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    grid % wall_dist(1:grid % n_cells),'TurbulentDistance')
  call Cgns_Mod_Write_Field_Seq(base, block, solution, field, grid, &
    grid % delta(1:grid % n_cells),'Cell_Delta')

  !----------------------!
  !   Pack in function   !
  !----------------------!

!!   go to base node
!      call cg_goto_f(index_file,index_base,ier,'end')
!!   write descriptor node (user can give any name)
!      text1='Supersonic vehicle with landing gear'
!      text2='M=4.6, Re=6 million'
!      textstring=text1//char(10)//text2
!      call cg_descriptor_write_f('Information',textstring,ier)

  ! Close DB
  call Cgns_Mod_Close_File_Seq

  print *, 'Successfully added fields to ', trim(problem_name)

  deallocate(cgns_base)

  ! Restore the name
  problem_name = store_name

  end subroutine
