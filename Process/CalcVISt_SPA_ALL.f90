!======================================================================!
  subroutine CalcVISt_SPA_ALL(n) 
!----------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                  !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, n 
  real    :: Xrat, Fv1, lf, Cs 
!======================================================================!

  if(SIMULA == DES_SPA) then
    do c = 1,NC
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      VISt(c)  = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  if(SIMULA == SPA_ALL) then
    do c = 1,NC
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      VISt(c)  = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  call Exchng(VISt)  

  RETURN

  end subroutine CalcVISt_SPA_ALL
