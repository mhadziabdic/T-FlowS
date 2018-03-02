!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   ~~~~~~~~~~~~~~                                                             !
!   For k-eps model :                                                          !
!                        2                                                     !
!   vis_t = Cmu * rho * k  * eps                                               ! 
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                           !
!            +            kappa                                                !
!   vis_tw = y  * vis_t ----------                                             ! 
!                       E * ln(y+)                                             !
!                                                                              !
!    For k-eps-v2f model :                                                     !
!                                                                              !
!    vis_t = CmuD * rho * Tsc  * vv                                            !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: ck, f_mu, re_t, y_star, pr_turb, ebf, pr_mol, beta                                        
!==============================================================================!

  !---------------------!
  !   k-epsilon model   !
  !---------------------!
  re_t    = 0.0
  f_mu    = 0.0 
  y_plus  = 0.0 
  pr_turb = 0.9

  ! Low-Re varaint
  if(turbulence_model_variant == LOW_RE) then
    do c = 1, grid % n_cells 
      re_t = kin % n(c)**2 / (viscosity * eps % n(c))
      f_mu = exp(-3.4/(1.0 + 0.02*re_t)**2) 
      vis_t(c) = f_mu * Cmu * density * kin % n(c)**2 / eps % n(c)
    end do

  ! High-Re varaint
  else
    do c = 1, grid % n_cells
      vis_t(c) = Cmu * density * kin % n(c)**2 / (eps % n(c) + TINY)
    end do
    if(ROUGH == NO) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
            ck = sqrt(tau_wall(c1))
            y_plus(c1) = density * ck * grid % wall_dist(c1) / viscosity 
            vis_wall(c1) = y_plus(c1) * viscosity *kappa  &
                         / log(Elog * y_plus(c1))
          end if
        end if
      end do
    else if(ROUGH == YES) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
            ck = sqrt(tau_wall(c1))
            y_plus(c1) = density*ck*(grid % wall_dist(c1)+Zo)/viscosity
            vis_wall(c1) = min(y_plus(c1) * viscosity * kappa                  &
                                      / log((grid % wall_dist(c1) + Zo) / Zo), &
                               1.0e+6 * viscosity)
          end if
        end if
      end do
    end if   
  end if

  ! Effective condctivity
  if(heat_transfer == YES) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          pr_mol = viscosity * capacity / conductivity
          con_wall(c1) = viscosity * capacity / pr_mol   
        end if
      end if
    end do
  end if
  
  !-----------------------!
  !   Hybrid PTIM model   !
  !-----------------------!
  if(turbulence_model == HYBRID_PITM) then
    do c = 1, grid % n_cells
      re_t = kin % n(c)*kin % n(c)/(viscosity*eps % n(c))

      y_star = (viscosity * eps % n(c))**0.25 * grid % wall_dist(c)/viscosity

      f_mu = (1.0 - exp(-y_star/14.0))**2.0*(1.0                              &
           + 5.0*exp(-(re_t/200.0)*(re_t/200.0))/re_t**0.75)


      f_mu = f_mu / ( 1.0 + exp(-y_star/5.0)**1.5/0.06 )
      f_mu = min(1.0,f_mu)

      vis_t(c) = f_mu * Cmu * density * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if

  call Comm_Mod_Exchange(grid, vis_t)  

  end subroutine
