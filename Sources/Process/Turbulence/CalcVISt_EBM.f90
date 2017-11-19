!==============================================================================!
  subroutine CalcVISt_EBM(grid) 
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
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: CK, Yplus, Fmu, Ret, yStar, Cmu_mod                                        
!==============================================================================!

  Ret   = 0.0
  Fmu   = 0.0 
  yPlus = 0.0 
  call CalcShear(U % n, V % n, W % n, Shear)

  do c = 1, grid % n_cells
    Kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

    Cmu_mod = max( -(uu % n(c) * Ux(c)                   &
                   + vv % n(c) * Vy(c)                   &
                   + ww % n(c) * Wz(c)                   &
                   + uv % n(c) * (Vx(c) + Uy(c))         &
                   + uw % n(c) * (Uz(c) + Wx(c))         &
                   + vw % n(c) * (Vz(c) + Wy(c)))        &
            / max(   Kin % n(c) * Kin % n(c)             &
                   / Eps%n(c)*Shear(c)*Shear(c), 1.0e-12), 0.0)

    Cmu_mod = min(0.12,Cmu_mod) 
    VISt(c) = Cmu_mod*DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
  end do

  call Exchange(grid, VISt)  

  end subroutine
