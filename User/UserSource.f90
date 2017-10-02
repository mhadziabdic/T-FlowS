!======================================================================!
  SUBROUTINE UserSource() 
!----------------------------------------------------------------------!
!   This subroutine extracts the heat in case when periodic boundaries !
!   are used in order to keep energy balans                            !
!   This source is part of the convection                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c
!----------------------------------------------------------------------!
! Description:                                                         !
! ~~~~~~~~~~~~                                                         !
!   Adds source to the energy equation.                                !
!======================================================================!
! 
!
!  /
! |
! |dT/dx * Ux * volume 
! |
!/
! dT/dx is derived from condition that there is no
! change of energy of the system. It means
! mass_flux * cp * dT - Q * Area = 0, Area = B*dx*Nwall,
!  dT/dx = Q*B*Nwall/(mass_flux*cp),
! where Q is heat flux through the wall, B is
! channel width and Nwall is number of heated walls.  
!
! For channel flow 
! b(c) = b(c) - U%n(c)*Q*Nwall*B/FLUXoX(material(c))  
!
! For pipe flow: dT/dz = D*pi*Q/FLUXoZ
!  
!  b(c)=b(c) - Qflux * W % n(c) / ( Wbulk(material(c)) + TINY ) * volume(c)
! 
!  In order to be consistent we did not use ideal d*pi but Area/L where
!  Area is total wall surface of pipe and L is lenght of pipe
!  AreaZ(1) is surface oposite to the flow stream
!  Qflux is calculated in CalBou
 

  if(JET == YES.or.BACKSTEP==YES.or.OTHER==YES) return
  
  if(CHANNEL==YES.and.PER_BC==YES) then
    do c=1,NC
      b(c) = b(c) -   1.0* Tflux * U % n(c) / (FLUXoX(material(c))) * volume(c)
!     b(c) = b(c) -   2.0*1.0* Tflux * U % n(c) / (FLUXoX(material(c))) * volume(c)
!     b(c)=b(c) -  2.0*3.14* Tflux * U % n(c) / (FLUXoX(material(c))) * volume(c)
    end do
  else if(PIPE == YES.and.PER_BC==YES) then
    do c=1,NC
      b(c)=b(c) -  2.0*3.1415926*0.005*W % n(c) / ( FLUXoZ(1) + TINY ) * volume(c)
    end do
  end if 

  return

  END SUBROUTINE UserSource
