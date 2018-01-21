!==============================================================================!
  subroutine Source_Vis_Spalart_Almaras(grid, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Computes the source terms in vis transport equation.                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Constants_Pro_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
  real    :: Xrat, Fv1, Fv2, Fw, SS, DistV, ProdV, R, GG, Dif
  real    :: dist
  real    :: phi_x(-grid % n_bnd_cells:grid % n_cells),  &
             phi_y(-grid % n_bnd_cells:grid % n_cells),  &
             phi_z(-grid % n_bnd_cells:grid % n_cells)
!==============================================================================!

  if(SIMULA == SPA_ALL) then
    do c = 1, grid % n_cells

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      Xrat  = VIS % n(c)/VISc
      Fv1   = Xrat**3/(Xrat**3 + Cvis1**3)
      Fv2   = 1.0 - Xrat/(1.0 + Xrat*Fv1)
      SS    = Vort(c) + VIS % n(c)*Fv2/(kappa**2*WallDs(c)**2)
      ProdV = Cb1 * DENc(material(c)) * SS * VIS % n(c)
      b(c)  = b(c) + ProdV * grid % vol(c)

      !----------------------------------!
      !   Compute the destruction term   !
      !----------------------------------!
      R     = VIS % n(c)/(SS * kappa**2 * WallDs(c)**2)
      GG    = R + Cw2*(R**6 - R)
      Fw    = GG*((1.0 + Cw3**6)/(GG**6 + Cw3**6))**(1.0/6.0)
      DistV = Cw1* DENc(material(c)) * Fw * (VIS % n(c)/WallDs(c)**2)
      A % val(A % dia(c)) = A % val(A % dia(c)) + DistV * grid % vol(c)
 
      !--------------------------------------------!
      !   Compute the first-order diffusion term   !
      !--------------------------------------------!
      Dif   = Cb2                                  &
            * DENc(material(c))                    &
            * (phi_x(c) + phi_y(c) + phi_z(c))**2  &
            / SIGMAv
      b(c)  = b(c) + Dif * grid % vol(c)
    end do

  else if(SIMULA == DES_SPA) then
    do c = 1, grid % n_cells

      ! What is 0.65 here?  A ghost number
      dist = min(WallDs(c),0.65 * grid % delta(c))

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      Xrat  = VIS % n(c)/VISc
      Fv1   = Xrat**3/(Xrat**3 + Cvis1**3)
      Fv2   = 1.0 - Xrat/(1.0 + Xrat*Fv1)
      SS    = Vort(c) + VIS % n(c)*Fv2/(kappa**2*dist**2)
      ProdV = Cb1 * DENc(material(c)) * SS * VIS % n(c)
      b(c)  = b(c) + ProdV * grid % vol(c)
      
      !-----------------------------------!
      !   Compute the destruction  term   !
      !-----------------------------------!
      R     = VIS % n(c)/(SS * kappa**2 * dist**2)
      GG    = R + Cw2*(R**6 - R)
      Fw    = GG*((1.0 + Cw3**6)/(GG**6 + Cw3**6))**(1.0/6.0)
      DistV = Cw1* DENc(material(c)) * Fw * (VIS % n(c)/dist**2)
      A % val(A % dia(c)) = A % val(A % dia(c)) + DistV * grid % vol(c)

      !--------------------------------------------!
      !   Compute the first-order diffusion term   !
      !--------------------------------------------!
      Dif   = Cb2                                  &
            * DENc(material(c))                    &
            * (phi_x(c) + phi_y(c) + phi_z(c))**2  &
            / SIGMAv
      b(c)  = b(c) + Dif * grid % vol(c)
    end do 
  end if

  end subroutine
