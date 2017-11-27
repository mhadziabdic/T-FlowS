!==============================================================================!
  subroutine Advection_Scheme(grid, &
                              phi_f, s,                           &
                              phi,                                &
                              phi_i, phi_j, phi_k, Di, Dj, Dk, &
                              blenda) 
!------------------------------------------------------------------------------!
!   Computes the value at the cell face using different convective  schemes.   !
!   In this subroutine I try to follow the nomenclature from Basara's and      !
!   Przulj's AIAA paper.                                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
  use Parameters_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi_f
  integer         :: s
  real            :: phi(-grid % n_bnd_cells:grid % n_cells)
  real            :: phi_i(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_j(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_k(-grid % n_bnd_cells:grid % n_cells)
  real            :: Di(grid % n_faces),  &
                     Dj(grid % n_faces),  &
                     Dk(grid % n_faces)
  integer         :: blenda
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, C, D
  real    :: fj ! flow oriented interpolation factor
  real    :: gD, gU, alfa, beta1, beta2 
  real    :: phij, phiU, phiUstar, rj, sign, GammaC, Beta
!==============================================================================!
!
!               Flux > 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  ==s=>  o    |    o    |----> xi
!   |    U    |    C    |    D    |         |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!
!               Flux < 0
!   +---------+---------+---------+---------+
!   |         |         |         |         |
!   |         |   c1    |   c2    |         |
!   |    o    |    o  <=s==  o    |    o    |----> xi
!   |         |    D    |    C    |    U    |
!   |         |         |         |         |
!   +---------+---------+---------+---------+   
!
!----------------------------------------------------------------------!

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  if(Flux(s) > 0.0) then ! goes from c1 to c2
    fj   = 1.0 - grid % f(s)
    C    = c1
    D    = c2
    sign = +1.0
  else ! Flux(s) < 0.0   ! goes from c2 to c1
    fj = grid % f(s)
    C    = c2
    D    = c1
    sign = -1.0
  end if

  if(Flux(s) > 0.0) then
    phiUstar = phi(D) - 2.0 * ( phi_i(C)*Di(s) &
                               +phi_j(C)*Dj(s) &
                               +phi_k(C)*Dk(s) )
  else
    phiUstar = phi(D) + 2.0 * ( phi_i(C)*Di(s) &
                               +phi_j(C)*Dj(s) &
                               +phi_k(C)*Dk(s) )
  end if

  phiU = max( phi_min(C), min(phiUstar, phi_max(C)) )

  rj = ( phi(C) - phiU ) / ( phi(D)-phi(C) + 1.0e-16 )

  gD = 0.5 * fj * (1.0+fj)
  gU = 0.5 * fj * (1.0-fj)

  if(blenda == CDS) then
    phij = fj
  else if(blenda == QUICK) then
    rj = ( phi(C) - phiU ) / ( phi(D)-phi(C) + 1.0e-12 )
    alfa = 0.0
    phij = (gD - alfa) + (gU + alfa) * rj
  else if(blenda == LUDS) then
    alfa = 0.5 * fj * (1+fj)
    phij = (gD - alfa) + (gU + alfa) * rj
  else if(blenda == MINMOD) then
    phij = fj * max(0.0, min(rj,1.0))
  else if(blenda == SMART) then
    beta1 = 3.0
    beta2 = 1.0
    phij = max( 0.0, min( (beta1-1.0)*rj, gD+gU*rj, beta2 ) )
  else if(blenda == AVL_SMART) then
    beta1 = 1.0 + fj*(2.0+fj) 
    beta2 = fj*(2.0-fj) 
    phij = max( 0.0, min( (beta1-1.0)*rj, gD+gU*rj, beta2 ) )
  else if(blenda == SUPERBEE) then
    phij = 0.5 * max( 0.0, min( 2.0*rj,1.0 ), min( rj,2.0 ) )
  else if(blenda == YES) then
    return
  end if

  phi_f = phi(C) + phij * sign * (phi(c2)-phi(c1))

  if(blenda == GAMMA) then
    Beta = 0.1

    if(Flux(s) > 0.0) then
      phiUstar = 1.0 - (phi(D) - phi(C))/(2.0 * ( phi_i(C)*Di(s) &
                                 +phi_j(C)*Dj(s) &
                                 +phi_k(C)*Dk(s)))
    else
      phiUstar = 1.0 + (phi(D) - phi(C))/(2.0 * ( phi_i(C)*Di(s) &
                                 +phi_j(C)*Dj(s) &
                                 +phi_k(C)*Dk(s)))
    end if

    GammaC = phiUstar/Beta

    if(phiUstar < Beta.and.phiUstar > 0.0) then
      phi_f = (1.0 - GammaC*(1.0 - grid % f(s))) * phi(C)   &
                   + GammaC*(1.0 - grid % f(s))  * phi(D)
    else if(phiUstar < 1.0.and.phiUstar >= Beta) then
       phi_f =        grid % f(s)  * phi(C)  &
             + (1.0 - grid % f(s)) * phi(D)
    else
      phi_f = phi(C)
    end if
  end if 

  end subroutine
