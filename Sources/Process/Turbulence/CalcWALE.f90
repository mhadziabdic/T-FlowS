!==============================================================================!
  subroutine CalcWALE(grid)
!------------------------------------------------------------------------------!
!  Compute SGS viscosity for LES by using WALE model.  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: ONE_THIRD
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
  real    :: SijdSijd(-grid % n_bnd_cells:grid % n_cells),  &
             She     (-grid % n_bnd_cells:grid % n_cells),  &
             Vor     (-grid % n_bnd_cells:grid % n_cells)
  integer :: c
  real    :: S11,  S22,  S33,  S12,  S13,  S23,  S21,  S31,  S32
  real    :: S11d, S22d, S33d, S12d, S13d, S23d, S21d, S31d, S32d
  real    :: V11,  V22,  V33,  V12,  V13,  V23,  V21,  V31,  V32
!==============================================================================!

  print *, '# I think there is a bug in this function (Bojan)'

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  do c = 1, grid % n_cells
    S11 = u % x(c)
    S22 = v % y(c)
    S33 = w % z(c)
    S12 = 0.5*(v % x(c) + u % y(c))
    S13 = 0.5*(u % z(c) + w % x(c))
    S23 = 0.5*(v % y(c) + w % y(c))
    S21 = S12
    S31 = S13
    S32 = S23

    V11 = 0.0   
    V22 = 0.0  
    V33 = 0.0  
    V12 = 0.5*(v % x(c) - u % y(c))
    V13 = 0.5*(u % z(c) - w % x(c))
    V23 = 0.5*(v % y(c) - w % y(c))
    V21 = -V12
    V31 = -V13
    V32 = -V23

    She(c) = 0.5 * shear(c) * shear(c)
    Vor(c) = 0.5 * vort(c) * vort(c)

    S11d =  S11*S11 + S12*S12 + S13*S13   &
         - (V11*V11 + V12*V12 + V13*V13)  &
         - ONE_THIRD * (She(c) - Vor(c))

    S22d =  S12*S12 + S22*S22 + S23*S23   &
         - (V12*V12 + V22*V22 + V23*V23)  &
         - ONE_THIRD * (She(c) - Vor(c))

    S33d =  S13*S13 + S23*S23 + S33*S33   &
         - (V13*V13 + V23*V23 + V33*V33)  &
         - ONE_THIRD * (She(c) - Vor(c))

    S12d = S11*S12 + S12*S22 + S13*S32 + (V11*V12 + V12*V22 + V13*V32) 
    S13d = S11*S13 + S12*S23 + S13*S33 + (V11*V13 + V12*V23 + V13*V33) 
    S23d = S21*S13 + S22*S23 + S23*S33 + (V21*V13 + V22*V23 + V23*V33) 

    S21d = S12d
    S31d = S13d
    S32d = S23d

    SijdSijd(c) = S11d*S11d + S22d*S22d + S33d*S33d  &
                + S12d*S12d + S13d*S13d + S23d*S23d
    
    WALEv(c) =  sqrt( abs (SijdSijd(c)**3) )          &
             / (sqrt( abs (She(c)     **5) ) +        &
                sqrt( sqrt(SijdSijd(c)**6) ) + TINY)
  end do 

  end subroutine
