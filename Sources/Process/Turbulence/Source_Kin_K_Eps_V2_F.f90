!==============================================================================!
  subroutine Source_Kin_K_Eps_V2_F(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in Kin transport equation.               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Parameters_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: Utan, UnorSq, Unor, UtotSq, dely, Stot
  real    :: lf, Gblend, Ustar, Ck, yPlus, Uplus, Pk_turb, Pk_vis
  real    :: EBF, EBF1
  real    :: ALPHA1, Lrans, Lsgs
!==============================================================================! 
!                                                                              !
!   Pk  is a production of turbulent kinematic energy.                         !
!                                                                              !
!   The form of Pk which is solving in this subroutine is :                    !
!                                                                              !
!   Pk = 2 * VISk { [ (dU/dx)**2 + (dV/dy)**2 + (dW/dz)**2 ] +                 ! 
!                   0.5( dV/dx + dU/dy )**2 +                                  !
!                   0.5( dW/dy + dV/dz )**2 +                                  !
!                   0.5( dU/dz + dW/dx )**2 }                                  !
!                                                                              !
!   Dimension of the Pk is [ kg m /s**3 ]                                      !
!   In kinetic energy equation exist two source terms which have form:         !
!                                                                              !
!     /           /                                                            !
!    |           |                                                             !
!    | Pk dV  -  | DENc Eps dV                                                 !
!    |           |                                                             !
!   /           /                                                              !
!                                                                              !
!------------------------------------------------------------------------------!

  if(SIMULA == HYB_ZETA) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      Lsgs  = 0.8*lf
      Lrans = 0.41*WallDs(c)
      ALPHA1 = max(1.0,Lrans/Lsgs)

      ! Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * grid % vol(c)

      ! Dissipation:
      if(ALPHA1 < 1.05) then
        A % val(A % dia(c)) = A % val(A % dia(c)) +          &
             DENc(material(c))*Eps % n(c)/(Kin%n(c) + TINY) * grid % vol(c)
      else
        A % val(A % dia(c)) = A % val(A % dia(c)) +                           &
                              DENc(material(c))   *                           &
                              min(ALPHA1**1.45 * Eps % n(c), Kin % n(c)**1.5  &
                            / (lf*0.01)) / (Kin % n(c) + TINY) * grid % vol(c)
      end if
      Pk(c) =  VISt(c) * Shear(c) * Shear(c)
    end do
  else
    do c = 1, grid % n_cells

      ! Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * grid % vol(c)
 
      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c)) +                             &
           DENc(material(c))*Eps % n(c)/(Kin % n(c)+TINY) * grid % vol(c)
      Pk(c) =  VISt(c) * Shear(c) * Shear(c) 
      if (BUOY == YES) then 
        buoyBeta(c) = 1.0
        Gbuoy(c) = -buoyBeta(c) * (grav_x * ut % n(c) +  &
                                   grav_y * vt % n(c) +  &
                                   grav_z * wt % n(c))
        b(c) = b(c) + max(0.0, Gbuoy(c) * grid % vol(c))
        A % val(A % dia(c)) = A % val(A % dia(c))                            &
                     + max(0.0,-Gbuoy(c)*grid % vol(c) / (Kin % n(c) + TINY))
      end if
    end do
  end if

  do s = 1, grid % n_faces
    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s)
    
    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
      if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

        ! Compute tangential velocity component
        UtotSq = U % n(c1) * U % n(c1) &
               + V % n(c1) * V % n(c1) &
               + W % n(c1) * W % n(c1)
        Unor = ( U % n(c1) * grid % sx(s)     &   
               + V % n(c1) * grid % sy(s)     &   
               + W % n(c1) * grid % sz(s) )   &   
               / sqrt(  grid % sx(s)*grid % sx(s)  &
                      + grid % sy(s)*grid % sy(s)  &
                      + grid % sz(s)*grid % sz(s))
        UnorSq = Unor*Unor

        if( UtotSq   >  UnorSq) then
          Utan = sqrt(UtotSq - UnorSq)
        else
          Utan = TINY
        end if
    
        if(Ynd(c1) > 3.0) then 
          if(ROUGH==NO) then
            TauWall(c1) = DENc(material(c1))*kappa*Uf(c1)*Utan    &   
                         /(log(Elog*Ynd(c1)))    
            Pk(c1) = TauWall(c1)*Uf(c1)/(kappa*WallDs(c1))
          else if(ROUGH==YES) then
            TauWall(c1) = DENc(material(c1))*kappa*Uf(c1)*Utan    &   
                           /(log((WallDs(c1)+Zo)/Zo))    
            Pk(c1) = TauWall(c1)*Uf(c1) / (kappa*(WallDs(c1)+Zo))
            Kin % n(c2) = TauWall(c1) / 0.09**0.5
          end if
          b(c1) = b(c1) + Pk(c1) * grid % vol(c1)
          b(c1) = b(c1) - VISt(c1) * Shear(c1) * Shear(c1) * grid % vol(c1)
        end if  
      end if  ! TypeBC(c2)==WALL or WALLFL
    end if    ! c2 < 0 
  end do

  end subroutine
