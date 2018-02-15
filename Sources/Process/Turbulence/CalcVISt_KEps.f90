!==============================================================================!
  subroutine Calcvist_Keps(grid) 
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
  use pro_mod
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
  real              :: CK, Yplus, Fmu, Ret, yStar, Prturb, EBF, Prmol, beta                                        
  character(len=80) :: heat_transfer
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
!==============================================================================!

  call Control_Mod_Heat_Transfer(heat_transfer)
  call Control_Mod_Turbulence_Model(turbulence_model)
  call Control_Mod_Turbulence_Model_Variant(turbulence_model_variant)

!======================================================================!
!  K-EPS model
!======================================================================!
  Ret   = 0.0
  Fmu   = 0.0 
  yPlus = 0.0 
  Prturb = 0.9

  if(turbulence_model_variant == 'HIGH_RE') then
    do c = 1, grid % n_cells
      vis_t(c) = Cmu * DENc(material(c)) * kin%n(c) * kin%n(c) / (eps % n(c)+1.0e-14)
    end do
    if(ROUGH==NO) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
          if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
            Ck = sqrt(TauWall(c1))
            yPlus = DENc(material(c1))*Ck*grid % wall_dist(c1)/VISc 
            VISwall(c1) = yPlus*VISc*kappa/LOG(Elog*yPlus)
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
            Ck = sqrt(TauWall(c1))
            yPlus = DENc(material(c1))*Ck*(grid % wall_dist(c1)+Zo)/VISc
            VISwall(c1) = min(yPlus*VISc*kappa/LOG((grid % wall_dist(c1)+Zo)/Zo),1.0e+6*VISc)
          end if
        end if
      end do
    end if   
  end if
  
  if(turbulence_model_variant == 'LOW_RE') then
    do c = 1, grid % n_cells 
      Ret = kin % n(c)*kin % n(c)/(VISc*eps % n(c))
      Fmu = exp(-3.4/(1.0 + 0.02*Ret)**2.0) 
      vis_t(c) = Fmu * Cmu * DENc(material(c)) * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if
  if(heat_transfer == 'YES') then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then  
        if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
          Prmol = VISc * CAPc(material(c1)) / CONc(material(c1))
          CONwall(c1) = VISc*CAPc(material(c1))/Prmol   
        end if
      end if
    end do
  end if   

  if(turbulence_model == 'HYB_PITM') then
    do c = 1, grid % n_cells
      Ret = kin % n(c)*kin % n(c)/(VISc*eps % n(c))

      yStar = (VISc * eps % n(c))**0.25 * grid % wall_dist(c)/VISc

      Fmu = (1.0 - exp(-yStar/14.0))**2.0*(1.0                              &
            + 5.0*exp(-(Ret/200.0)*(Ret/200.0))/Ret**0.75)


      Fmu = Fmu / ( 1.0 + exp(-yStar/5.0)**1.5/0.06 )
      Fmu = min(1.0,Fmu)

      vis_t(c) = Fmu * Cmu * DENc(material(c)) * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if

  call Exchange(grid, vis_t)  

  end subroutine
