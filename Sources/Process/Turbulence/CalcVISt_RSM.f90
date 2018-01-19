!==============================================================================!
  subroutine CalcVISt_RSM(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models (EBM and HJ).              !
!   If hybrid option is used turbulent diffusivity is modeled by VISt.         !
!   Otherwise, VISt is used as false diffusion in order to increase            !
!   stability of computation.                                                  !
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
  real    :: Cmu_mod                                        
!==============================================================================!

  call Compute_Shear_And_Vorticity(grid)

  if(SIMULA == HJ) then
    do c=1, grid % n_cells
      Kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

      Cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(Kin % n(c)**2 / Eps_tot(c) * Shear(c)**2, 1.0e-12), 0.0)

      Cmu_mod = min(0.12,Cmu_mod) 
      VISt(c) = Cmu_mod * DENc(material(c)) * Kin % n(c)**2 / Eps_tot(c)
    end do 
  else if(SIMULA==EBM) then
    do c=1, grid % n_cells
      Kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

      Cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(Kin % n(c)**2 / Eps % n(c) * Shear(c)**2, 1.0e-12), 0.0)

      Cmu_mod = min(0.12,Cmu_mod)
      VISt(c) = Cmu_mod*DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
    end do
  end if

  call Exchange(grid, VISt)  

  end subroutine
