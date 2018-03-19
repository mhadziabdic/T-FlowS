!==============================================================================!
  subroutine Source_Kin_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation.               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
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
  real              :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real              :: lf
  real              :: alpha1, l_rans, l_sgs
!==============================================================================! 
!                                                                              !
!   p_kin  is a production of turbulent kinematic energy.                      !
!                                                                              !
!   The form of p_kin which is solving in this subroutine is :                 !
!                                                                              !
!   p_kin = 2 * VISk { [ (dU/dx)**2 + (dV/dy)**2 + (dW/dz)**2 ] +              ! 
!                   0.5( dV/dx + dU/dy )**2 +                                  !
!                   0.5( dW/dy + dV/dz )**2 +                                  !
!                   0.5( dU/dz + dW/dx )**2 }                                  !
!                                                                              !
!   Dimension of the p_kin is [ kg m /s**3 ]                                   !
!   In kinetic energy equation exist two source terms which have form:         !
!                                                                              !
!     /           /                                                            !
!    |           |                                                             !
!    | p_kin dV  -  | DENc eps dV                                              !
!    |           |                                                             !
!   /           /                                                              !
!                                                                              !
!------------------------------------------------------------------------------!

  if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      l_sgs  = 0.8*lf
      l_rans = 0.41*grid % wall_dist(c)
      alpha1 = max(1.0,l_rans/l_sgs)

      ! Production:
      b(c) = b(c) + vis_t(c) * shear(c) * shear(c) * grid % vol(c)

      ! Dissipation:
      if(alpha1 < 1.05) then
        A % val(A % dia(c)) = A % val(A % dia(c)) + &
             density * eps % n(c)/(kin%n(c) + TINY) * grid % vol(c)
      else
        A % val(A % dia(c)) = A % val(A % dia(c)) +                           &
                              density *                                       &
                              min(alpha1**1.45 * eps % n(c), kin % n(c)**1.5  &
                            / (lf*0.01)) / (kin % n(c) + TINY) * grid % vol(c)
      end if
      p_kin(c) =  vis_t(c) * shear(c) * shear(c)
    end do
  else
    do c = 1, grid % n_cells

      ! Production:
      p_kin(c) =  vis_t(c)/density * shear(c) * shear(c)
      b(c) = b(c) + density * p_kin(c) * grid % vol(c)
 
      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c)) + &
           density * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)

      if (buoyancy == YES) then 
        buoyBeta(c) = 1.0
        Gbuoy(c) = -buoyBeta(c) * (grav_x * ut % n(c) +  &
                                   grav_y * vt % n(c) +  &
                                   grav_z * wt % n(c)) * density
        b(c) = b(c) + max(0.0, Gbuoy(c) * grid % vol(c))
        A % val(A % dia(c)) = A % val(A % dia(c))  &
                     + max(0.0,-Gbuoy(c)*grid % vol(c) / (kin % n(c) + TINY))
      end if
    end do
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then

        ! Compute tangential velocity component
        u_tot_sq = u % n(c1) * u % n(c1) &
                 + v % n(c1) * v % n(c1) &
                 + w % n(c1) * w % n(c1)
        u_nor = ( u % n(c1) * grid % sx(s)     &
                + v % n(c1) * grid % sy(s)     &
                + w % n(c1) * grid % sz(s) )   &
                / sqrt(  grid % sx(s)*grid % sx(s)  &
                       + grid % sy(s)*grid % sy(s)  &
                       + grid % sz(s)*grid % sz(s))
        u_nor_sq = u_nor*u_nor

        if( u_tot_sq   >  u_nor_sq) then
          u_tan = sqrt(u_tot_sq - u_nor_sq)
        else
          u_tan = TINY
        end if
    
        if(y_plus(c1) > 3.0) then 
          if (ROUGH.eq.NO) then
            ! dimesions here are incorrect !!!
            tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                         /(log(Elog*y_plus(c1)))
            p_kin(c1) = tau_wall(c1)*u_tau(c1)/ &
              (kappa*grid % wall_dist(c1))/density
          else if (ROUGH.eq.YES) then
            tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                           /(log((grid % wall_dist(c1)+Zo)/Zo))
            p_kin(c1) = tau_wall(c1)*u_tau(c1) / &
              (kappa*(grid % wall_dist(c1)+Zo))/density
            kin % n(c2) = tau_wall(c1) / 0.09**0.5
          end if
          b(c1) = b(c1) + density * p_kin(c1) * grid % vol(c1)
          b(c1) = b(c1) - vis_t(c1) * shear(c1) * shear(c1) * grid % vol(c1)
        end if  
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL or WALLFL
    end if    ! c2 < 0 
  end do

  end subroutine
