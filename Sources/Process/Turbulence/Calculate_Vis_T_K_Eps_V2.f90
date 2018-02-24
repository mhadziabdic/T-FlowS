!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps_V2(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use allp_mod
  use Flow_Mod
  use Comm_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: u_nor_sq, u_nor, u_tot_sq, Cmu1, beta, pr_mol, pr_turb
  real    :: lf, g_blend, u_star, c_k, y_plus, u_plus, ebf
!==============================================================================!

  call Time_And_Length_Scale(grid)

  if(turbulence_model == K_EPS_V2) then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v2%n(c)*Tsc(c)
    end do
  else if(turbulence_model == K_EPS_ZETA_F) then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v2%n(c)*kin % n(c)*Tsc(c)
    end do
  else if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v2%n(c)*kin % n(c)*Tsc(c)
      vis_t_eff(c) = max(vis_t(c),vis_t_sgs(c))
    end do
    call Comm_Mod_Exchange(grid, vis_t_eff)  
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then

        Uf(c1)  = Cmu**0.25*kin%n(c1)**0.5

        if(ROUGH == YES) then
          Ynd(c1) = (grid % wall_dist(c1)+Zo)*Uf(c1)/viscosity
        else if(ROUGH == NO) then
          Ynd(c1) = grid % wall_dist(c1)*Uf(c1)/viscosity 
        end if

        Uf(c1)  = Cmu**0.25*kin%n(c1)**0.5
        Ynd(c1) = grid % wall_dist(c1)*Uf(c1)/viscosity 
        g_blend  = 0.01*Ynd(c1)**4.0/(1.0+5.0*Ynd(c1))

        y_plus = max(Ynd(c1),0.13)  
        u_plus = log(y_plus*Elog)/(kappa)
   
        if(y_plus< 3.0) then
          VISwall(c1) = vis_t(c1) + viscosity  
        else
          VISwall(c1) = Ynd(c1)*viscosity/(y_plus*exp(-1.0*g_blend) &
                        +u_plus*exp(-1.0/g_blend) + TINY)
        end if

        if(ROUGH == YES) then
          u_plus = log((grid % wall_dist(c1)+Zo)/Zo)/(kappa + TINY) + TINY
          VISwall(c1) = min(y_plus*viscosity*kappa/LOG((grid % wall_dist(c1)+Zo)/Zo),1.0e+6*viscosity)
       end if

        if(heat_transfer == YES) then
          pr_turb = 1.0 / ( 0.5882 + 0.228*vis_t(c1)/viscosity   &
                 - 0.0441 * (vis_t(c1)/viscosity)**2.0  &
                 * (1.0 - exp(-5.165*viscosity/(vis_t(c1)+tiny))) )
          pr_mol = viscosity * capacity / conductivity
          beta = 9.24 * ((pr_mol/pr_turb)**0.75 - 1.0) * (1.0 + 0.28 * exp(-0.007*pr_mol/pr_turb))
          ebf = 0.01 * (pr_mol*y_plus)**4.0 / (1.0 + 5.0 * pr_mol**3.0 * y_plus) + TINY
          CONwall(c1) = y_plus*viscosity*capacity/(y_plus*pr_mol*exp(-1.0 * ebf) &
                      + (u_plus + beta)*pr_turb*exp(-1.0/ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL or WALLFL
    end if    ! c2 < 0
  end do
 
  call Comm_Mod_Exchange(grid, vis_t)  
  call Comm_Mod_Exchange(grid, VISwall)  

  end subroutine
