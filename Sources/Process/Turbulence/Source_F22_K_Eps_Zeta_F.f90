!==============================================================================!
  subroutine Source_F22_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation  equation                     !
!   for vi2 and imposing  boundary condition for f22                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, j
  real              :: sor_11, f22hg
  real              :: A0
!==============================================================================!
!                                                                              !
!  The form of source terms are :                                              !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22hg*dV ; f22hg - f22hg homogenious is placed in a source              !
!    |                     coefficients b(c)                                   !
!   /                                                                          !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/kin(c) - 2.0/3.0)/Tsc(c)     &             !
!              + 2.0*Cv2*p_kin(c)/(3.0*kin(c))                                !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22*dV ;this term is placed in a diagonal of coefficient matrix         !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!                                                                              !
!  Dimensions of certain variables                                             !
!                                                                              !
!     Tsc            [s]                                                       !
!     kin            [m^2/s^2]                                                 !
!     eps            [m^3/s^2]                                                 !
!     vi2            [m^2/s^2]                                                 !
!     f22            [-]                                                       !
!     Lsc            [m]                                                       !
!                                                                              !
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

 ! Source term f22hg
 if(turbulence_model == K_EPS_ZETA_F .or.  &
    turbulence_model == HYBRID_K_EPS_ZETA_F) then 
   do c = 1, grid % n_cells
     f22hg = (1.0 - Cf_1 - 0.65 * p_kin(c)  &
           / (eps  % n(c) + TINY))          &
           * (zeta % n(c) - TWO_THIRDS)     &
           / (Tsc(c) + TINY)                &
           + 0.0085 * p_kin(c) / (kin % n(c) + TINY)
     b(c) = b(c) + f22hg * grid % vol(c) / (Lsc(c)**2 + TINY) 
   end do

   ! Source term f22hg
   do c = 1, grid % n_cells
     sor_11 = grid % vol(c)/(Lsc(c)**2 + TINY)
     A % val(A % dia(c)) = A % val(A % dia(c)) + sor_11 
   end do

   ! Imposing boundary condition for f22 on the wall
   do s = 1, grid % n_faces
     c1=grid % faces_c(1,s)
     c2=grid % faces_c(2,s)
     if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER ) then
       if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
          Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then

         f22 % n(c2) = -2.0 * viscosity * zeta % n(c1)     &
                     / grid % wall_dist(c1)**2

        ! Fill in a source coefficients

        ! Linearization of the near wall terms helps to get more  
         ! stable solution, especially for small wall distance.
         A0 = f_coef(s)
         b(c1) = b(c1) + A0 * f22 % n(c2)
       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do
 end if 

 end subroutine