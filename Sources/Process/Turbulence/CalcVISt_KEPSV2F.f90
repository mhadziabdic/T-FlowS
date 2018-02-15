!==============================================================================!
  subroutine Calcvist_KepsV2F(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
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
  real              :: UnorSq, Unor, UtotSq, Cmu1, beta, Prmol, Prturb
  real              :: lf, Gblend, Ustar, Ck, yPlus, Uplus, EBF
  character(len=80) :: heat_transfer
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
!==============================================================================!

  call Time_And_Length_Scale(grid)

  call Control_Mod_Heat_Transfer(heat_transfer)
  call Control_Mod_Turbulence_Model(turbulence_model)

  if(turbulence_model == 'K_EPS_VV') then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v_2%n(c)*Tsc(c)
    end do
  else if(turbulence_model == 'ZETA') then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v_2%n(c)*kin % n(c)*Tsc(c)
    end do
  else if(turbulence_model == 'HYB_ZETA') then
    do c = 1, grid % n_cells
      vis_t(c) = CmuD*v_2%n(c)*kin % n(c)*Tsc(c)
      vis_t_eff(c) = max(vis_t(c),vis_t_sgs(c))
    end do
    call Exchange(grid, vis_t_eff)  
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then

        Uf(c1)  = Cmu**0.25*kin%n(c1)**0.5

        if(ROUGH == YES) then
          Ynd(c1) = (grid % wall_dist(c1)+Zo)*Uf(c1)/VISc
        else if(ROUGH == NO) then
          Ynd(c1) = grid % wall_dist(c1)*Uf(c1)/VISc 
        end if

        Uf(c1)  = Cmu**0.25*kin%n(c1)**0.5
        Ynd(c1) = grid % wall_dist(c1)*Uf(c1)/VISc 
        Gblend  = 0.01*Ynd(c1)**4.0/(1.0+5.0*Ynd(c1))

        Yplus = max(Ynd(c1),0.13)  
        Uplus = log(Yplus*Elog)/(kappa)
   
        if(Yplus< 3.0) then
          VISwall(c1) = vis_t(c1) + VISc  
        else
          VISwall(c1) = Ynd(c1)*VISc/(Yplus*exp(-1.0*Gblend) &
                        +Uplus*exp(-1.0/Gblend) + TINY)
        end if

        if(ROUGH == YES) then
          Uplus = log((grid % wall_dist(c1)+Zo)/Zo)/(kappa + TINY) + TINY
          VISwall(c1) = min(Yplus*VISc*kappa/LOG((grid % wall_dist(c1)+Zo)/Zo),1.0e+6*VISc)
       end if

        if(heat_transfer == 'YES') then
          Prturb = 1.0 / ( 0.5882 + 0.228*vis_t(c1)/VISc   &
                 - 0.0441 * (vis_t(c1)/VISc)**2.0  &
                 * (1.0 - exp(-5.165*VISc/(vis_t(c1)+tiny))) )
          Prmol = VISc * CAPc(material(c1)) / CONc(material(c1))
          beta = 9.24 * ((Prmol/Prturb)**0.75 - 1.0) * (1.0 + 0.28 * exp(-0.007*Prmol/Prturb))
          EBF = 0.01 * (Prmol*Yplus)**4.0 / (1.0 + 5.0 * Prmol**3.0 * Yplus) + TINY
          CONwall(c1) = Yplus*VISc*CAPc(material(c1))/(Yplus*Prmol*exp(-1.0 * EBF) &
                      + (Uplus + beta)*Prturb*exp(-1.0/EBF) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL or WALLFL
    end if    ! c2 < 0
  end do
 
  call Exchange(grid, vis_t)  
  call Exchange(grid, VISwall)  

  end subroutine
