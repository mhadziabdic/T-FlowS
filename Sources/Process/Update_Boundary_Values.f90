!==============================================================================!
  subroutine Update_Boundary_Values(grid)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, s
  real              :: qx, qy, qz, Nx, Ny, Nz, Stot, CONeff, EBF, Yplus
  real              :: Prmol, beta, Uplus, Prt
  character(len=80) :: heat_transfer
  character(len=80) :: turbulence_model
!==============================================================================!

  Area  = 0.0
  Tflux = 0.0

  call Control_Mod_Heat_Transfer(heat_transfer)
  call Control_Mod_Turbulence_Model(turbulence_model)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    !------------------------------------------------!
    !   Outflow (and inflow, if needed) boundaries   !
    !------------------------------------------------!

    ! On the boundary perform the extrapolation
    if(c2  < 0) then

      ! Extrapolate velocities on the outflow boundary 
      ! SYMMETRY is intentionally not treated here because I wanted to
      ! be sure that is handled only via graPHI and NewUVW functions)
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  &
          Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE ) then
        U % n(c2) = U % n(c1)
        V % n(c2) = V % n(c1)
        W % n(c2) = W % n(c1)
        if(heat_transfer == 'YES') t % n(c2) = t % n(c1)
      end if

      if( Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY ) then
        U % n(c2) = U % n(c1)
        V % n(c2) = V % n(c1)
        W % n(c2) = W % n(c1)
        if(heat_transfer == 'YES') t % n(c2) = t % n(c1)
      end if

      ! Spalart Allmaras
      if(turbulence_model == 'SPA_ALL' .or.  &
         turbulence_model == 'DES_SPA') then
        if ( Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  & 
             Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE ) then
          VIS % n(c2) = VIS % n(c1) 
        end if
      end if

      if(turbulence_model == 'EBM' .or.  &
         turbulence_model == 'HJ') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          uu % n(c2) = 0.0 
          vv % n(c2) = 0.0 
          ww % n(c2) = 0.0 
          uv % n(c2) = 0.0 
          uw % n(c2) = 0.0 
          vw % n(c2) = 0.0 
          KIN % n(c2) = 0.0 
          if(turbulence_model == 'EBM') f22% n(c2) = 0.0 
        end if
      end if

      ! k-epsilon-v^2
      if(turbulence_model == 'K_EPS_VV' .or.  &
         turbulence_model == 'ZETA'     .or.  &
         turbulence_model == 'HYB_ZETA') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          v_2 % n(c2) = v_2 % n(c1)
          f22 % n(c2) = f22 % n(c1)
        end if

        !  if (Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW) then
        !    f22 % n(c2) = f22 % n(c1)
        !  end if
      end if 

      ! k-epsilon
      if(turbulence_model == 'K_EPS') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
        end if
      end if 

      if(turbulence_model == 'EBM' .or.  &
         turbulence_model == 'HJ') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          uu % n(c2) = uu % n(c1)
          vv % n(c2) = vv % n(c1)
          ww % n(c2) = ww % n(c1)
          uv % n(c2) = uv % n(c1)
          uw % n(c2) = uw % n(c1)
          vw % n(c2) = vw % n(c1)
          if(turbulence_model == 'EBM') f22 % n(c2) = f22 % n(c1)
        end if
      end if 

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(heat_transfer == 'YES') then
        Prt = 0.9
        if(turbulence_model/='LES' .or.  &
           turbulence_model/='DNS') then
          Prt = 1.0                           &
              / ( 0.5882 + 0.228*(vis_t(c1)    &
              / (VISc+1.0e-12)) - 0.0441      &
              * (vis_t(c1)/(VISc+1.0e-12))**2  &
              * (1.0 - exp(-5.165*( VISc/(vis_t(c1)+1.0e-12) ))))
        end if
        Stot = sqrt(  grid % sx(s)*grid % sx(s)  &
                    + grid % sy(s)*grid % sy(s)  &
                    + grid % sz(s)*grid % sz(s))
        Nx = grid % sx(s)/Stot
        Ny = grid % sy(s)/Stot
        Nz = grid % sz(s)/Stot
        qx = t % q(c2) * Nx 
        qy = t % q(c2) * Ny
        qz = t % q(c2) * Nz
        CONeff = CONc(material(c1))                 &
               + CAPc(material(c1))*vis_t(c1)/Prt
        if(turbulence_model == 'ZETA' .or.  &
           turbulence_model == 'K_EPS') then
          Yplus = max(Cmu25 * sqrt(kin%n(c1)) * grid % wall_dist(c1)/VISc,0.12)
          Uplus = log(Yplus*Elog) / (kappa + TINY) + TINY
          Prmol = VISc / CONc(material(c1))
          beta = 9.24 * ((Prmol/Prt)**0.75 - 1.0)  &
                     * (1.0 + 0.28 * exp(-0.007*Prmol/Prt))
          EBF = 0.01 * (Prmol*Yplus)**4.0          &
                     / (1.0 + 5.0 * Prmol**3 * Yplus) + TINY
          CONwall(c1) = Yplus * VISc * CAPc(material(c1))   &
                      / (Yplus * Prmol * exp(-1.0 * EBF)    &
                      + (Uplus + beta) * Prt * exp(-1.0 / EBF) + TINY)
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
            t% n(c2) = t % n(c1) + Prt / CAPc(material(c1))  &
                     * (  qx * grid % dx(s)                  &
                        + qy * grid % dy(s)                  &
                        + qz * grid % dz(s))                 &
                     / CONwall(c1)
            Tflux = t % q(c2)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * CONeff  &
                      / (  Nx * grid % dx(s)                &
                         + Ny * grid % dy(s)                &
                         + Nz * grid % dz(s) )
            Tflux = t % q(c2)
          end if
        else
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
            t% n(c2) = t % n(c1) + Prt / (CONeff + TINY)   &
                    * (  qx * grid % dx(s)                 &
                       + qy * grid % dy(s)                 &
                       + qz * grid % dz(s) ) 
            Tflux = t % q(c2) 
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * CONeff  &
                      / (  Nx * grid % dx(s)                &
                         + Ny * grid % dy(s)                &
                         + Nz * grid % dz(s) )
            Tflux = t % q(c2) 
          end if
        end if
      end if

      !---------------------!
      !   Copy boundaries   !
      !---------------------!
      if(grid % bnd_cond % copy_c(c2) /= 0) then
        U % n(c2) = U % n(grid % bnd_cond % copy_c(c2))
        V % n(c2) = V % n(grid % bnd_cond % copy_c(c2))
        W % n(c2) = W % n(grid % bnd_cond % copy_c(c2))

        if(heat_transfer == 'YES')  &
          t % n(c2) = t % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model == 'SPA_ALL' .or.   &
           turbulence_model == 'DES_SPA') &
          VIS % n(c2) = VIS % n(grid % bnd_cond % copy_c(c2)) 

        if(turbulence_model == 'K_EPS_VV' .or.  &
           turbulence_model == 'ZETA'     .or.  &
           turbulence_model == 'HYB_ZETA') then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          v_2 % n(c2) = v_2 % n(grid % bnd_cond % copy_c(c2))
          f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if ! 'K_EPS_VV'

        if(turbulence_model == 'K_EPS') then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
        end if ! 'K_EPS'

        if(turbulence_model == 'EBM' .or.  &
           turbulence_model == 'HJ') then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          uu % n(c2)  = uu % n(grid % bnd_cond % copy_c(c2))
          vv % n(c2)  = vv % n(grid % bnd_cond % copy_c(c2))
          ww % n(c2)  = ww % n(grid % bnd_cond % copy_c(c2))
          uv % n(c2)  = uv % n(grid % bnd_cond % copy_c(c2))
          uw % n(c2)  = uw % n(grid % bnd_cond % copy_c(c2))
          vw % n(c2)  = vw % n(grid % bnd_cond % copy_c(c2))
          f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if ! 'EBM' 
      end if
    end if
  end do

  !  if(heat_transfer == 'YES') then 
  !    call GloSum(Tflux)
  !    call GloSum(Area)
  !    call wait
  !   
  !    Tflux  = Tflux/(Area+TINY)
  !    Qflux  = Tflux * Area/(4.0*AreaZ(1))
  !  end if 
  !  print *, Tflux
  !  This is done for pipe flow
  !  In order to be consistent we did not use ideal d*pi but Area/L where
  !  Area is total wall surface of pipe and L is lenght of pipe
  !  AreaZ(1) is surface in direction of to the flow stream

  end subroutine
