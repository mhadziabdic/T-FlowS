!==============================================================================!
  subroutine Source_Eps_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Calculates source terms in equation of dissipation of turbulent energy     !
!   and imposes boundary condition                                             ! 
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
  integer :: c, s, c1, c2,j 
  real    :: Esor, Ce_11, EBF
  real    :: Utan, UnorSq, Unor, UtotSq, EpsWall, EpsHom
  real    :: p_kin_turb, p_kin_vis, y_pl
!==============================================================================!
!   In dissipation of turbulent kinetic energy equation exist two              !
!   source terms which have form:                                              !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | ((Cv_e1 * p_kinE - Cv_11 DENc * eps) / Tsc) * dV                           !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!   First, positive , source term is solved and added to source  coefficient   !
!   b(c) on right hand side.  Second, negative, source term is added to main   !
!   diagonal left hand side coefficient matrix in order to increase stability  !
!   of solver.  It is nessesary to calculate coefficient Cv_11 using kin,      !
!   Cv_e2, vi2 and coefficient A1                                              !
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    do c = 1, grid % n_cells 
      Esor = grid % vol(c)/(Tsc(c)+TINY)
      Ce_11 = Ce1*(1.0 + alpha * ( 1.0/(zeta % n(c)+TINY) ))    
      b(c) = b(c) + Ce_11*p_kin(c)*Esor
 
      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + Ce2*Esor*density
    end do                   
  end if

  ! Imposing a boundary condition on wall for eps 

  do s = 1, grid % n_faces
    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s)
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER ) then
     if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
        Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
        UtotSq = U % n(c1) * U % n(c1) &
               + V % n(c1) * V % n(c1) &
               + W % n(c1) * W % n(c1) 
        Unor = ( U % n(c1) * grid % sx(s)     &
               + V % n(c1) * grid % sy(s)     &
               + W % n(c1) * grid % sz(s) )   &
               / sqrt(  grid % sx(s)*grid % sx(s)  &
                      + grid % sy(s)*grid % sy(s)  &
                      + grid % sz(s)*grid % sz(s) )  
        UnorSq = Unor*Unor

        if( UtotSq   >  UnorSq) then
          Utan = sqrt(UtotSq - UnorSq)
        else
          Utan = TINY
        end if

        EpsWall = 2.0 * viscosity/density * kin % n(c1) / grid % wall_dist(c1)**2
        EpsHom = Cmu75 * kin % n(c1)**1.5 / (grid % wall_dist(c1) * kappa)
        u_tau(c1) = Cmu25 * kin % n(c1)**0.5


        p_kin_turb = Cmu75 * kin % n(c1)**1.5 / (grid % wall_dist(c1) * kappa) * density
        p_kin_vis = vis_t(c1)*Shear(c1)*Shear(c1) * density ! standard

        ! Kader
        ! EBF = 0.01 * y_pl**4.0 / (1.0 + 5.0*y_pl) + TINY           !original
        ! Pro = p_kin_vis * exp(-1.0 * EBF) + p_kin_turb * exp(-1.0 / EBF) 
        ! BL_EPS = min((p_kin_turb * exp(-1.0 / EBF))/Pro,1.0)

        y_pl = Cmu25 * sqrt(kin % n(c1)) * grid % wall_dist(c1) / (viscosity/density + TINY) !standard
        EBF  = 0.001*y_pl**4.0/(1.0+y_pl)
        eps % n(c1) = EpsWall * exp(-1.0 * EBF) + EpsHom * exp(-1.0 / EBF)
        
        if(ROUGH == YES) then
          eps % n(c1) = Cmu75 * kin % n(c1)**1.5 / ((grid % wall_dist(c1)) * kappa)
        end if

        ! Adjusting coefficient to fix eps value in near wall calls
        do j = A % row(c1), A % row(c1+1) - 1
          A % val(j) = 0.0
        end do
        b(c1) = eps % n(c1) * density
        A % val(A % dia(c1)) = 1.0 * density
     end if
    end if
  end do  

!@if(turbulence_model == K_EPS_V2_F) then
!@  do c = 1, grid % n_cells
!@    Esor = grid % vol(c)/Tsc(c)
!@    Ce_11 = Ce1 * (1.0 + alpha * sqrt(kin % n(c) / (v2 % n(c)) + TINY))
!@    b(c) = b(c) + Ce_11 * p_kin(c) * Esor
!@
!@    ! Fill in a diagonal of coefficient matrix
!@    A % val(A % dia(c)) =  A % val(A % dia(c)) + Ce2 * Esor * density
!@  end do
!@end if 

  end subroutine
