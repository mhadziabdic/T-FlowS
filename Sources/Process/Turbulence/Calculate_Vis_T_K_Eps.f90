!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   ~~~~~~~~~~~~~~                                                             !
!   For k-eps model :                                                          !
!                        2                                                     !
!   vis_t = Cmu * rho * K  * eps                                               ! 
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                           !
!            +          kappa                                                  !
!   vis_tw = y  * vis_t ----------                                             ! 
!                     E * ln(y+)                                               !
!                                                                              !
!    For k-eps-v2f model :                                                     !
!                                                                              !
!    vis_t = CmuD * rho * Tsc  * vv                                            !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, c1, c2, s
  real              :: CK, Yplus, Fmu, Ret, yStar, Prturb, EBF, pr_mol, beta                                        
!==============================================================================!

  !---------------------!
  !   k-epsilon model   !
  !---------------------!
  Ret   = 0.0
  Fmu   = 0.0 
  yPlus = 0.0 
  Prturb = 0.9

  ! High-Re varaint
  if(turbulence_model_variant == HIGH_RE) then
    do c = 1, grid % n_cells
      vis_t(c) = Cmu * density * kin%n(c) * kin%n(c) / (eps % n(c)+1.0e-14)
    end do
    if(ROUGH==NO) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
          if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
            Ck = sqrt(tau_wall(c1))
            yPlus = density*Ck*grid % wall_dist(c1)/viscosity 
            VISwall(c1) = yPlus*viscosity*kappa/log(Elog*yPlus)
          end if
        end if
      end do
    else if(ROUGH==YES) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
            Ck = sqrt(tau_wall(c1))
            yPlus = density*Ck*(grid % wall_dist(c1)+Zo)/viscosity
            VISwall(c1) = min(yPlus*viscosity*kappa/log((grid % wall_dist(c1)+Zo)/Zo),1.0e+6*viscosity)
          end if
        end if
      end do
    end if   
  end if
  
  ! Low-Re varaint
  if(turbulence_model_variant == LOW_RE) then
    do c = 1, grid % n_cells 
      Ret = kin % n(c)*kin % n(c)/(viscosity*eps % n(c))
      Fmu = exp(-3.4/(1.0 + 0.02*Ret)**2.0) 
      vis_t(c) = Fmu * Cmu * density * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if
  if(heat_transfer == YES) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          pr_mol = viscosity * capacity / conductivity
          CONwall(c1) = viscosity * capacity / pr_mol   
        end if
      end if
    end do
  end if   

  !-----------------------!
  !   Hybrid PTIM model   !
  !-----------------------!
  if(turbulence_model == HYBRID_PITM) then
    do c = 1, grid % n_cells
      Ret = kin % n(c)*kin % n(c)/(viscosity*eps % n(c))

      yStar = (viscosity * eps % n(c))**0.25 * grid % wall_dist(c)/viscosity

      Fmu = (1.0 - exp(-yStar/14.0))**2.0*(1.0                              &
            + 5.0*exp(-(Ret/200.0)*(Ret/200.0))/Ret**0.75)


      Fmu = Fmu / ( 1.0 + exp(-yStar/5.0)**1.5/0.06 )
      Fmu = min(1.0,Fmu)

      vis_t(c) = Fmu * Cmu * density * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if

  call Exchange(grid, vis_t)  

  end subroutine
