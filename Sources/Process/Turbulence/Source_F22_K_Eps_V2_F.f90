!==============================================================================!
  subroutine Source_F22_K_Eps_V2_F(grid)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation  equation                     !
!   for vi2 and imposing  boundary condition for f22                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer s, c, c1, c2, j
  real    Sor11,  f22hg
  real    A0
!------------------------------------------------------------------------------!
!                                                                              !
!  The form of source terms are :                                              !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22hg*dV ; f22hg - f22hg homogenious is placed in a source              !
!    |                     coefficients b(c)                                   !
!   /                                                                          !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/Kin(c) - 2.0/3.0)/Tsc(c)     &             !
!              + 2.0*Cv_2*Pk(c)/(3.0*Kin(c))                                   !
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
!     Kin            [m^2/s^2]                                                 !
!     Eps            [m^3/s^2]                                                 !
!     vi2            [m^2/s^2]                                                 !
!     f22            [-]                                                       !
!     Lsc            [m]                                                       !
!==============================================================================!

  call Scale()

 ! Source term f22hg
 if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then 
   do c = 1,NC
     f22hg = (1.0 - Cf_1 - 0.65*Pk(c)   &
           / (Eps % n(c) + TINY))       &
           * (v_2 % n(c) - 2.0 / 3.0)   &
           / (Tsc(c) + TINY)            &
           + 0.0085 * Pk(c) / (Kin % n(c) + TINY)
     b(c) = b(c) + f22hg * volume(c) / (Lsc(c)**2 + TINY) 
   end do

   ! Source term f22hg
   do c = 1, NC
     Sor11 = volume(c)/(Lsc(c)**2 + TINY)
     A % val(A % dia(c)) = A % val(A % dia(c)) + Sor11 
   end do

   ! Imposing boundary condition for f22 on the wall
   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then


         f22 % n(c2) = -2.0 * VISc * v_2 % n(c1)     &
                     / WallDs(c1)**2

        ! Fill in a source coefficients

        ! Linearization of the near wall terms helps to get more  
         ! stable solution, especially for small wall distance.
         A0 = Scoef(s) 
         b(c1) = b(c1) + A0*f22%n(c2)
       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do
 else if(SIMULA == K_EPS_VV) then
   do c = 1,NC
     f22hg = (1.0 - Cf_1)*(v_2 % n(c)/(Kin % n(c)+TINY) - 2.0/3.0)/(Tsc(c)+TINY)   &
                + Cf_2*Pk(c)/(Kin % n(c)+TINY)
     b(c) = b(c) + f22hg*volume(c)/(Lsc(c)+TINY)**2
     Sor11 = volume(c)/Lsc(c)**2
     A % val(A % dia(c)) = A % val(A % dia(c)) + Sor11
   end do

   ! Imposing boundary condition for f22 on the wall
   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
         f22 % n(c2) = -20.0*VISc**2*v_2 % n(c1)/                     &
                       (WallDs(c1)**4*Eps % n(c2))

        ! Fill in a source coefficients

        ! Linearization of the near wall terms helps to get more  
         ! stable solution, especially for small wall distance.
         A0 = Scoef(s)
         b(c1) = b(c1) + A0*f22%n(c2)
       end if   ! end if of BC=wall
     end if    ! end if of c2<0
   end do
 end if 

 end subroutine
