!==============================================================================!
  subroutine Update_Boundary_Values(grid)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, s
  real              :: qx, qy, qz, nx, ny, nz, Stot, CONeff, ebf, y_plus
  real              :: Prmol, beta, u_plus, Prt
!==============================================================================!

  area  = 0.0
  Tflux = 0.0

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
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        if(heat_transfer == YES) t % n(c2) = t % n(c1)
      end if

      if( Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY ) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        if(heat_transfer == YES) t % n(c2) = t % n(c1)
      end if

      ! Spalart Allmaras
      if(turbulence_model == SPALART_ALLMARAS .or.  &
         turbulence_model == DES_SPALART) then
        if ( Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  & 
             Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE ) then
          vis % n(c2) = vis % n(c1) 
        end if
      end if

      if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
         turbulence_model == HANJALIC_JAKIRLIC) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          uu % n(c2) = 0.0 
          vv % n(c2) = 0.0 
          ww % n(c2) = 0.0 
          uv % n(c2) = 0.0 
          uw % n(c2) = 0.0 
          vw % n(c2) = 0.0 
          KIN % n(c2) = 0.0 
          if(turbulence_model == REYNOLDS_STRESS_MODEL) f22% n(c2) = 0.0 
        end if
      end if

      ! k-epsilon-v^2
      if(turbulence_model == K_EPS_V2 .or.  &
         turbulence_model == K_EPS_ZETA_F     .or.  &
         turbulence_model == HYBRID_K_EPS_ZETA_F) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          v2  % n(c2) = v2  % n(c1)
          f22 % n(c2) = f22 % n(c1)
        end if

        !  if (Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW) then
        !    f22 % n(c2) = f22 % n(c1)
        !  end if
      end if 

      ! k-epsilon
      if(turbulence_model == K_EPS) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
        end if
      end if 

      if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
         turbulence_model == HANJALIC_JAKIRLIC) then
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
          if(turbulence_model == REYNOLDS_STRESS_MODEL) f22 % n(c2) = f22 % n(c1)
        end if
      end if 

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(heat_transfer == YES) then
        Prt = 0.9
        if(turbulence_model /= LES .or.  &
           turbulence_model /= DNS) then
          Prt = 1.0                           &
              / ( 0.5882 + 0.228*(vis_t(c1)    &
              / (viscosity+1.0e-12)) - 0.0441      &
              * (vis_t(c1)/(viscosity+1.0e-12))**2  &
              * (1.0 - exp(-5.165*( viscosity/(vis_t(c1)+1.0e-12) ))))
        end if
        Stot = sqrt(  grid % sx(s)*grid % sx(s)  &
                    + grid % sy(s)*grid % sy(s)  &
                    + grid % sz(s)*grid % sz(s))
        nx = grid % sx(s)/Stot
        ny = grid % sy(s)/Stot
        nz = grid % sz(s)/Stot
        qx = t % q(c2) * nx 
        qy = t % q(c2) * ny
        qz = t % q(c2) * nz
        CONeff = conductivity                 &
               + capacity*vis_t(c1)/Prt
        if(turbulence_model == K_EPS_ZETA_F .or.  &
           turbulence_model == K_EPS) then
          y_plus = max(Cmu25 * sqrt(kin%n(c1)) * grid % wall_dist(c1)/viscosity,0.12)
          u_plus = log(y_plus*Elog) / (kappa + TINY) + TINY
          Prmol = viscosity / conductivity
          beta = 9.24 * ((Prmol/Prt)**0.75 - 1.0)  &
                     * (1.0 + 0.28 * exp(-0.007*Prmol/Prt))
          ebf = 0.01 * (Prmol*y_plus)**4.0          &
                     / (1.0 + 5.0 * Prmol**3 * y_plus) + TINY
          con_wall(c1) = y_plus * viscosity * capacity   &
                      / (y_plus * Prmol * exp(-1.0 * ebf)    &
                      + (u_plus + beta) * Prt * exp(-1.0 / ebf) + TINY)
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
            t% n(c2) = t % n(c1) + Prt / capacity  &
                     * (  qx * grid % dx(s)                  &
                        + qy * grid % dy(s)                  &
                        + qz * grid % dz(s))                 &
                     / con_wall(c1)
            Tflux = t % q(c2)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * CONeff  &
                      / (  nx * grid % dx(s)                &
                         + ny * grid % dy(s)                &
                         + nz * grid % dz(s) )
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
                      / (  nx * grid % dx(s)                &
                         + ny * grid % dy(s)                &
                         + nz * grid % dz(s) )
            Tflux = t % q(c2) 
          end if
        end if
      end if

      !---------------------!
      !   Copy boundaries   !
      !---------------------!
      if(grid % bnd_cond % copy_c(c2) /= 0) then
        u % n(c2) = u % n(grid % bnd_cond % copy_c(c2))
        v % n(c2) = v % n(grid % bnd_cond % copy_c(c2))
        w % n(c2) = w % n(grid % bnd_cond % copy_c(c2))

        if(heat_transfer == YES)  &
          t % n(c2) = t % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model == SPALART_ALLMARAS .or.   &
           turbulence_model == DES_SPALART) &
          vis % n(c2) = vis % n(grid % bnd_cond % copy_c(c2)) 

        if(turbulence_model == K_EPS_V2 .or.  &
           turbulence_model == K_EPS_ZETA_F     .or.  &
           turbulence_model == HYBRID_K_EPS_ZETA_F) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          v2  % n(c2) = v2  % n(grid % bnd_cond % copy_c(c2))
          f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if 

        if(turbulence_model == K_EPS) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
        end if 

        if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
           turbulence_model == HANJALIC_JAKIRLIC) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          uu % n(c2)  = uu  % n(grid % bnd_cond % copy_c(c2))
          vv % n(c2)  = vv  % n(grid % bnd_cond % copy_c(c2))
          ww % n(c2)  = ww  % n(grid % bnd_cond % copy_c(c2))
          uv % n(c2)  = uv  % n(grid % bnd_cond % copy_c(c2))
          uw % n(c2)  = uw  % n(grid % bnd_cond % copy_c(c2))
          vw % n(c2)  = vw  % n(grid % bnd_cond % copy_c(c2))
          f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if 
      end if
    end if
  end do

  end subroutine
