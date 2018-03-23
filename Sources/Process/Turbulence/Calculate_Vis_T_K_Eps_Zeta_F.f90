!=============================================================================!
  subroutine Calculate_Vis_T_K_Eps_Zeta_F(grid) 
!-----------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.               !
!---------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
  use Work_Mod,    only: kin_visc_t => r_cell_01
  use Work_Mod,    only: visc_ratio => r_cell_02
!-----------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: beta, pr_mol, pr_turb
  real    :: g_blend, y_pl, u_plus, ebf
  real    :: kin_visc
!=============================================================================!

  call Time_And_Length_Scale(grid)

  if(turbulence_model == K_EPS_ZETA_F) then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD * zeta % n(c) * kin % n(c) * Tsc(c) * density
    end do
  else if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD * zeta % n(c) * kin % n(c) * Tsc(c)
      vis_t_eff(c) = max(vis_t(c), vis_t_sgs(c))
    end do
    call Comm_Mod_Exchange(grid, vis_t_eff)  
  end if

  ! kinematic viscosities
  kin_visc   = viscosity/density
  kin_visc_t = vis_t    /density
  visc_ratio = vis_t    /viscosity

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tau(c1)  = Cmu**0.25*kin%n(c1)**0.5

        if(ROUGH == YES) then
          y_plus(c1) = (grid % wall_dist(c1)+Zo)*u_tau(c1)/kin_visc
        else if(ROUGH == NO) then
          y_plus(c1) = grid % wall_dist(c1)*u_tau(c1)/kin_visc 
        end if

        g_blend  = 0.01*y_plus(c1)**4.0/(1.0+5.0*y_plus(c1))

        y_pl = max(y_plus(c1),0.13)  
        u_plus = log(y_pl*Elog)/kappa
   
        if(y_pl < 3.0) then
          ! kg / (m * s) !
          vis_wall(c1) = vis_t(c1) + viscosity
        else
          vis_wall(c1) = y_plus(c1)*viscosity/ &
            (y_pl*exp(-1.0*g_blend) + u_plus*exp(-1.0/g_blend) + TINY)
        end if

        if(ROUGH == YES) then
          u_plus = log((grid % wall_dist(c1)+Zo)/Zo)/(kappa + TINY) + TINY
          vis_wall(c1) = min(y_pl*viscosity*kappa/ &
            LOG((grid % wall_dist(c1)+Zo)/Zo),1.0e+6*kin_visc)
        end if

        if(heat_transfer == YES) then
          pr_turb = 1.0 / ( 0.5882 + 0.228*visc_ratio(c1) &
                  - 0.0441 * visc_ratio(c1)**2.0          &
                  * (1.0 - exp(-5.165/(visc_ratio(c1)+TINY)) ) )

          pr_mol = kin_visc * capacity / conductivity
          beta = 9.24 * ((pr_mol/pr_turb)**0.75 - 1.0) * &
            (1.0 + 0.28 * exp(-0.007*pr_mol/pr_turb))
          ebf = 0.01 * (pr_mol*y_pl)**4.0 / &
            (1.0 + 5.0 * pr_mol**3.0 * y_pl) + TINY
          con_wall(c1) = y_pl*kin_visc*capacity/(y_pl*pr_mol* &
            exp(-1.0 * ebf) + (u_plus + beta)*pr_turb*exp(-1.0/ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do
 
  call Comm_Mod_Exchange(grid, vis_t)  
  call Comm_Mod_Exchange(grid, vis_wall)  

  end subroutine