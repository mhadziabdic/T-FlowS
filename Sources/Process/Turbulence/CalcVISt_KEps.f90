!==============================================================================!
  subroutine CalcVISt_Keps(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   ~~~~~~~~~~~~~~                                                             !
!   For k-eps model :                                                          !
!                       2                                                      !
!   VISt = Cmu * rho * K  * Eps                                                ! 
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                           !
!            +          kappa                                                  !
!   VIStw = y  * VISt ----------                                               ! 
!                     E * ln(y+)                                               !
!                                                                              !
!    For k-eps-v2f model :                                                     !
!                                                                              !
!    VISt = CmuD * rho * Tsc  * vv                                             !
!                                                                              !
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
  real    :: CK, Yplus, Fmu, Ret, yStar, Prturb, EBF, Prmol, beta                                        
!==============================================================================!

!======================================================================!
!  K-EPS model
!======================================================================!
  Ret   = 0.0
  Fmu   = 0.0 
  yPlus = 0.0 
  Prturb = 0.9

  if(MODE == HRe) then
    do c = 1, grid % n_cells
      VISt(c) = Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / (Eps % n(c)+1.0e-14)
    end do
    if(ROUGH==NO) then
      do s = 1, grid % n_faces
        c1 = SideC(1,s)
        c2 = SideC(2,s)
        if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then  
          if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
            Ck = sqrt(TauWall(c1))
            yPlus = DENc(material(c1))*Ck*WallDs(c1)/VISc 
            VISwall(c1) = yPlus*VISc*kappa/LOG(Elog*yPlus)
          end if
        end if
      end do
    else if(ROUGH==YES) then
      do s = 1, grid % n_faces
        c1 = SideC(1,s)
        c2 = SideC(2,s)
        if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
          if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
            Ck = sqrt(TauWall(c1))
            yPlus = DENc(material(c1))*Ck*(WallDs(c1)+Zo)/VISc
            VISwall(c1) = min(yPlus*VISc*kappa/LOG((WallDs(c1)+Zo)/Zo),1.0e+6*VISc)
          end if
        end if
      end do
    end if   
  end if
  
  if(MODE==LRe) then
    do c = 1, grid % n_cells 
      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))
      Fmu = exp(-3.4/(1.0 + 0.02*Ret)**2.0) 
      VISt(c) = Fmu * Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
    end do
  end if
  if(HOT == YES) then
    do s = 1, grid % n_faces
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then  
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          Prmol = VISc * CAPc(material(c1)) / CONc(material(c1))
          CONwall(c1) = VISc*CAPc(material(c1))/Prmol   
        end if
      end if
    end do
  end if   

  if(SIMULA==HYB_PITM) then
    do c = 1, grid % n_cells
      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))

      yStar = (VISc * Eps % n(c))**0.25 * WallDs(c)/VISc

      Fmu = (1.0 - exp(-yStar/14.0))**2.0*(1.0                              &
            + 5.0*exp(-(Ret/200.0)*(Ret/200.0))/Ret**0.75)


      Fmu = Fmu / ( 1.0 + exp(-yStar/5.0)**1.5/0.06 )
      Fmu = min(1.0,Fmu)

      VISt(c) = Fmu * Cmu * DENc(material(c)) * Kin%n(c) * Kin%n(c) / Eps % n(c)
    end do
  end if

  call Exchange(grid, VISt)  

  end subroutine
