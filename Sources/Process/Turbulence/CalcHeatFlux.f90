!======================================================================!
  SUBROUTINE CalcHeatFlux(grid)
!----------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                     !
!                                                                      !
!   Authors: Muhamed Hadziabdic                                        ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Parameters_Mod
  use Work_Mod, only: T_x => r_cell_01,  &
                      T_y => r_cell_02,  &
                      T_z => r_cell_03    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
  real    :: beta, Pr
!==============================================================================!

  call GraPhi(grid, T % n, 1, T_x, .TRUE.)
  call GraPhi(grid, T % n, 2, T_y, .TRUE.)
  call GraPhi(grid, T % n, 3, T_z, .TRUE.)

!-------------------------------------------!
!    Compute the sources in the interior    !
!-------------------------------------------!
  Pr  = 0.71
  beta = 1.0
  Prt = 0.9

  if(SIMULA==K_EPS.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA.or.&
     SIMULA==DES_SPA) then
    do c = 1, grid % n_cells
      Prt = 1.0/( 0.5882 + 0.228*(VISt(c)/(VISc+1.0e-12)) - 0.0441 * &   
            (VISt(c)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*&
           ( VISc/(VISt(c)+1.0e-12) ))) )
      ut%n(c) = -VISt(c)/Prt * T_x(c)
      vt%n(c) = -VISt(c)/Prt * T_y(c)
      wt%n(c) = -VISt(c)/Prt * T_z(c)
     
      if(BUOY == YES) then 
        ut%n(c) = min(0.01*Tref,ut%n(c))
        ut%n(c) = max(-0.01*Tref,ut%n(c))
        vt%n(c) = min(0.01*Tref,vt%n(c))
        vt%n(c) = max(-0.01*Tref,vt%n(c))
        wt%n(c) = min(0.01*Tref,wt%n(c))
        wt%n(c) = max(-0.01*Tref,wt%n(c))
        Pbuoy(c) = -beta*(grav_x*ut%n(c) + grav_y*vt%n(c) + grav_z*wt%n(c))
        Pbuoy(c) = max(Pbuoy(c),0.0)
      end if
    end do
  else if(SIMULA==EBM.or.SIMULA==HJ) then
    do c = 1, grid % n_cells
      Prt = 1.0/( 0.5882 + 0.228*(VISt(c)/(VISc+1.0e-12)) - 0.0441 * &   
            (VISt(c)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*&
           ( VISc/(VISt(c)+1.0e-12) ))) )
      ut%n(c) =  -0.22*Tsc(c) * (uu % n(c) * T_x(c) + &
                 uv % n(c) * T_y(c) + uw % n(c) * T_z(c))
      vt%n(c) =  -0.22*Tsc(c) * (uv % n(c) * T_x(c) + &
                 vv % n(c) * T_y(c) + vw % n(c) * T_z(c))
      wt%n(c) =  -0.22*Tsc(c) * (uw % n(c) * T_x(c) + &
                 vw % n(c) * T_y(c) + ww % n(c) * T_z(c))

      if(BUOY == YES) then 
        ut%n(c) = min(0.01*Tref,ut%n(c))
        ut%n(c) = max(-0.01*Tref,ut%n(c))
        vt%n(c) = min(0.01*Tref,vt%n(c))
        vt%n(c) = max(-0.01*Tref,vt%n(c))
        wt%n(c) = min(0.01*Tref,wt%n(c))
        wt%n(c) = max(-0.01*Tref,wt%n(c))
        Pbuoy(c) = -beta*(grav_x*ut%n(c) + grav_y*vt%n(c) + grav_z*wt%n(c))
        Pbuoy(c) = max(Pbuoy(c),0.0)
      end if
    end do
  end if
  RETURN

  END SUBROUTINE CalcHeatFlux
