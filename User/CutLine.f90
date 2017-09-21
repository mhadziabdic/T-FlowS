!======================================================================!
  SUBROUTINE CutLine(namAut)
!----------------------------------------------------------------------!
!  Purpose: Writes values of variables in cutlines aligned with coord. !
!           axes for arbitrary mesh.                                   !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER, OPTIONAL :: namAut*(*)
!------------------------------[Calling]-------------------------------!
  REAL      :: Dist
!-------------------------------[Locals]-------------------------------!
  CHARACTER :: answer*140,  ReadLine1*140, XYZ*5
  LOGICAL   :: Hscale, Here(NC)
  INTEGER   :: s, c1, c2, Ncl, Cln, cp, j, l, Nwc, UC, DC, dir
  INTEGER   :: Xdir, Ydir, Zdir, RealM, BulkM, UtauM, PlusM, CLnorm
  REAL      :: Jcut, Kcut, Weight
  REAL      :: Ic(-NbC:NC), Jc(-NbC:NC), Kc(-NbC:NC)
  REAL      :: u_tau, u_tau1, u_tau2, Re_tau, Re_h, BulkVel, Yplus, TauWall
  REAL      :: Unorm, Hnorm, Tnorm, viscos, Flog, Qflux, Tinf, Twall
  REAL      :: Icut(NC),  Ucl(NC),  Vcl(NC),  Wcl(NC),  Tcl(NC),       &
               Kin_cl(NC), Eps_cl(NC),  V_2_cl(NC),  F22_cl(NC),       &
               Pk_cl(NC),  UiUj_cl(NC), UiUk_cl(NC), VISt_cl(NC)
!--------------------------------[CVS]---------------------------------!
!  '$Id: CutLine.f90,v 1.2 2017/08/31 22:42:35 mhadziabdic Exp $'/
!  '$Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/User/CutLine.f90,v $'/
!======================================================================!


!-----Parameter numbers------------------------------------------------!
  Xdir  = 1
  Ydir  = 2
  Zdir  = 3
  RealM = 4
  BulkM = 5
  UtauM = 6
  PlusM = 7

!-----Reading number of cutlines---------------------------------------!

!--during calculation
  if(PRESENT(namAut)) then

!--name
    namCut = namAut

!--go on in T-Rex.cmn file
    call ReadC(7,inp,tn,ts,te)
    ReadLine1(1:len_trim(inp)) = inp(1:len_trim(inp))
    call ReadC(7,inp,tn,ts,te)
    call ReadC(7,inp,tn,ts,te)
    call ReadC(7,inp,tn,ts,te)

!--number of CL
    read(inp(ts(2):te(2)),*) Ncl

!--from T-Rex.cmn file
  else
    if(This < 2) write(*,*) '# Writing T-REX cutline file name [skip cancels]:'

!--name
    call ReadC(7,inp,tn,ts,te)
    read(inp(ts(1):te(1)),'(A80)') namCut
    answer = namCut
    call ToUppr(answer)
    if(answer == 'SKIP') then
      if(This < 2) write(*,*) '# Writing T-REX cutline skiped!!'
      RETURN
    end if

!--number of CL
    read(inp(ts(2):te(2)),*) Ncl
  end if
  namCut( len_trim(namCut)+1:len_trim(namCut)+8 ) = "-cut.000"


!=====From here everything is in one pass for each CL==================!
  do Cln=1,Ncl

!--reading parameters
    Hscale = .FALSE.
    call ReadC(7,inp,tn,ts,te)

!--direction
    read(inp(ts(1):te(1)),'(A8)') answer
    call ToUppr(answer)
    if(answer == 'X') then
      dir = Xdir
    else if(answer == 'Y') then
      dir = Ydir
    else if(answer == 'Z') then
      dir = Zdir
    else
      write(*,*) '@CutLine: Unknown direction! No CL created!!'
      RETURN
    end if

!--constant coordinates
    read(inp(ts(2):te(2)),*) Jcut
    read(inp(ts(3):te(3)),*) Kcut

!--normalization mode
    read(inp(ts(4):te(4)),'(A8)') answer
    call ToUppr(answer)
    if(answer == 'REAL') then 
      CLnorm = RealM
      if(tn == 5) then
        read(inp(ts(5):te(5)),*) Hnorm
        Hscale = .TRUE.
      end if
    else if(answer == 'UTAU') then
      CLnorm = UtauM
      if(tn == 5) then
        read(inp(ts(5):te(5)),*) Hnorm
        Hscale = .TRUE.
      end if
    else if(answer == 'BULK') then
      CLnorm = BulkM
      if(tn == 5) then
        read(inp(ts(5):te(5)),*) Hnorm
        Hscale = .TRUE.
      end if
    else if(answer == 'ADD') then
      CLnorm = PlusM
      if(tn == 5) then
        read(inp(ts(5):te(5)),*) Hnorm
        Hscale = .TRUE.
      end if
    else
      write(*,*) '@CutLine: Wrong normalization mode! No CL created!!'
      RETURN
    end if

!--REAL normalization for DNS
    if(SIMULA==DNS) then
      CLnorm = RealM
      write(*,*) '@CutLine: DNS - Normalization mode set to REAL!'
    end if

!--opening file
    l = len_trim(namCut)
    if (Cln < 10) then
      write(namcut(l:l),'(I1)') Cln
    else if (Cln < 100) then
      write(namcut(l-1:l),'(I2)') Cln
    else
      write(namcut(l-2:l),'(I3)') Cln
    end if
    if(This < 2) write(*,*) '# NOW CREATING FILE: ', namCut
    open (100+Cln, FILE=namCut)

!--write file header
    if(dir == Xdir) XYZ = 'x'
    if(dir == Ydir) XYZ = 'y'
    if(dir == Zdir) XYZ = 'z'

    write(100+Cln,'(A10,I3,A2,A1,A14,F7.3,F7.3,A2)')            &
    ' # Cutline', Cln,'; ', XYZ, '-axis (const: ', Jcut, Kcut,')'

    if(CLnorm /= PlusM) then
      if(CLnorm==BulkM) write(100+Cln,'(A)') ' # Non-normalized'
      if(CLnorm==UtauM) write(100+Cln,'(A)') ' # U_tau normalized'
      if(CLnorm==BulkM) write(100+Cln,'(A)') ' # Bulk normalized'
      write (100+Cln,*) '#'
      write (100+Cln,'(A3,A1,A39)') ' # ', XYZ,                          &
      ', U, V, W, T, Kin, Eps, Pk, uv, v_2, f22'
    else
      write(100+Cln,'(A)') ' # Specified normalization'
      write (100+Cln,*) '#'
      write (100+Cln,'(A3,A1,A23)') ' # ', XYZ, ', Cf, Y+, U/Ub, V/Ub, St'
    end if


!-----Find CL cells and their variable values--------------------------!


!--define appropriate direction
    if (dir == Xdir) then
      Ic = xc
      Jc = yc
      Kc = zc
      Qflux = -1.0 * HUGE
    else if (dir == Ydir) then
      Ic = yc
      Jc = xc
      Kc = zc
      Qflux = -1.0 * HUGE
    else if (dir == Zdir) then
      Ic = zc
      Jc = xc
      Kc = yc
      Qflux = -1.0 * HUGE
    end if

!--initialization
    Here = .FALSE.

!--browse through faces
    do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!--HeatFlux
      if(c2 < 0 .and. TypeBC(c2)==WALLFL) then
        if(HOT==YES) Qflux = max( Qflux, T%q(c2) )
      end if

!--Jcut in between c1 and c2
      if( c2>0 .and. ( (Jc(c1) >= Jcut .and. Jc(c2) <= Jcut) .or.   &
                       (Jc(c1) <= Jcut .and. Jc(c2) >= Jcut) ) ) then

        if( abs(Jc(c1)-Jcut) < abs(Jc(c2)-Jcut) ) then
          do j=Acol(c1),Acol(c1+1)-1
            if(Arow(j) > 0 .and. Arow(j) /= c1) then
!--Kcut in between c1 and Arow
              if( (Kc(c1) >= Kcut .and. Kc(Arow(j)) <= Kcut) .or. &
                  (Kc(c1) <= Kcut .and. Kc(Arow(j)) >= Kcut) ) then
                Here(c1) = .TRUE.
                UC = c1
                DC = c2
              end if
            end if
          end do
        else
          do j=Acol(c2),Acol(c2+1)-1
            if(Arow(j) > 0 .and. Arow(j) /= c2) then
!--Kcut in between c2 and Arow
              if( (Kc(c2) >= Kcut .and. Kc(Arow(j)) <= Kcut) .or. &
                  (Kc(c2) <= Kcut .and. Kc(Arow(j)) >= Kcut) ) then
                Here(c2) = .TRUE.
                UC = c2
                DC = c1
              end if
            end if !end if Arow
          end do   !end do j=
        end if     !end if c1 or c2 closer to Jcut
      end if       !end if Jc1>Jcut>Jc2

      if( IsNearWall(c1) ) then
        Nwc = c1
      else if( IsNearWall(c2) ) then
        Nwc = c2
      end if

      if( Here(c1) .or. Here(c2) ) then
        Weight = Dist( Ic(DC), Jc(DC), Kc(DC), Ic(DC), Jcut,   Kcut ) &
               / Dist( Ic(UC), Jc(UC), Kc(UC), Ic(DC), Jc(DC), Kc(DC) )

        Icut(UC) = Ic(UC) * Weight + Ic(DC) * (1.0 - Weight)

        Ucl(UC) = U % n(UC) * Weight + U % n(DC) * (1.0 - Weight)

        Vcl(UC) = V % n(UC) * Weight + V % n(DC) * (1.0 - Weight)

        Wcl(UC) = W % n(UC) * Weight + W % n(DC) * (1.0 - Weight)

        if(HOT==YES)                                            &
        Tcl(UC) = T % n(UC) * Weight + T % n(DC) * (1.0 - Weight)

        if(SIMULA==K_EPS .or. SIMULA==K_EPS_VV) then
          Kin_cl(UC) = Kin%n(UC) * Weight + Kin%n(DC) * (1.0 - Weight)

          Eps_cl(UC) = Eps%n(UC) * Weight + Eps%n(DC) * (1.0 - Weight)

          if(SIMULA==K_EPS_VV) then
            V_2_cl(UC) = V_2%n(UC) * Weight + V_2%n(DC) * (1.0 - Weight)

            F22_cl(UC) = F22%n(UC) * Weight + F22%n(DC) * (1.0 - Weight)
          end if

          VISt_cl(UC) = VISt(UC) * Weight + VISt(DC) * (1.0 - Weight)

          Pk_cl(UC) = Pk(UC) * Weight + Pk(DC) * (1.0 - Weight)

          if(dir == Xdir) then
            UiUj_cl(UC) = (Uy(UC) + Vx(UC)) * Weight       &
                        + (Uy(DC) + Vx(DC)) * (1.0 - Weight)

            UiUk_cl(UC) = (Uz(UC) + Wx(UC)) * Weight       &
                        + (Uz(DC) + Wx(DC)) * (1.0 - Weight)

          else if(dir == Ydir) then
            UiUj_cl(UC) = (Uy(UC) + Vx(UC)) * Weight       &
                        + (Uy(DC) + Vx(DC)) * (1.0 - Weight)

            UiUk_cl(UC) = (Vz(UC) + Wy(UC)) * Weight       &
                        + (Vz(DC) + Wy(DC)) * (1.0 - Weight)

          else if(dir == Zdir) then
            UiUj_cl(UC) = (Uz(UC) + Wx(UC)) * Weight       &
                        + (Uz(DC) + Wx(DC)) * (1.0 - Weight)

            UiUk_cl(UC) = (Vz(UC) + Wy(UC)) * Weight       &
                        + (Vz(DC) + Wy(DC)) * (1.0 - Weight)
          end if

          UiUj_cl(UC) = UiUj_cl(UC) * VISt_cl(UC)
          UiUk_cl(UC) = UiUk_cl(UC) * VISt_cl(UC)

          if(SIMULA==K_EPS .and. UC==Nwc) then
            if(WFType==ANL) then
              UiUj_cl(Nwc) = VISwall(Nwc)*Utan(Nwc)/WallDs(Nwc)
              UiUk_cl(Nwc) = VISwall(Nwc)*Utan(Nwc)/WallDs(Nwc)
            else
              UiUj_cl(Nwc) = VISt_cl(Nwc)*(u_tau/(kappa*WallDs(Nwc)))
              UiUj_cl(Nwc) = VISt_cl(Nwc)*(u_tau/(kappa*WallDs(Nwc)))
            end if
          end if

        end if    !end if SIMULA==K_EPS
      end if      !end if(Here c1,c2)
    end do        !end do s=1,NS

!-----Defining normalizating values------------------------------------!

!--calculate u_tau
    if(SIMULA==K_EPS) then
      u_tau1 = Cmu25 * sqrt( Kin_cl(Nwc) )
      u_tau2 = abs( Utan(Nwc)*kappa/log(Elog*WallDs(Nwc)*u_tau1/VISc) )
      u_tau = u_tau1
!      u_tau = u_tau2
!      u_tau = sqrt( u_tau1 * u_tau2 )
!      u_tau = sqrt(abs( VISwall(Nwc)*Utan(Nwc)/WallDs(Nwc) ))
!      u_tau = sqrt(abs(PdropX(material(1))))
    end if
    if(SIMULA==K_EPS_VV) then
      if(MODE==WF) then
!       u_tau = Ufmean(Nwc)
        if(ZETA==YES) then                                           !Zeta!--New u_tau
          u_tau = sqrt(sqrt(CmuD * Kin_cl(Nwc)**2.0 * v_2_cl(Nwc)))  !Zeta!--New u_tau
        else
          u_tau = sqrt(sqrt(CmuD * Kin_cl(Nwc) * v_2_cl(Nwc)))
        end if
      else
        u_tau = sqrt(abs( VISc*Utan(Nwc)/WallDs(Nwc) ))
      end if
    end if

!--normalizating values
    viscos =  1.0
    Unorm  =  1.0
!    Tnorm  = -1.0
    Tnorm  = -20.0
    if( .NOT. Hscale ) Hnorm = 1.0
    BulkVel = max( FluxX(1)/AreaX(1), &
                   FluxY(1)/AreaY(1), &
                   FluxZ(1)/AreaZ(1)  )

    if(CLnorm==UtauM) then       !U_tau normalization
      if(u_tau==0.0) then
        if(This < 2) write(*,*) '@CutLine: Problem with U_tau!!'
      else
        viscos= VISc
        Unorm = u_tau
        if(HOT==YES) Tnorm = 1.0 / (DENc(material(Nwc))*CAPc(material(Nwc))*u_tau)
        if( .NOT. Hscale ) Hnorm = VISc / u_tau
      end if
    else if(CLnorm==BulkM) then       !BulkVel normalization
      viscos= VISc
      Unorm = BulkVel
    else if(CLnorm==PlusM) then       !ADD - additional normalization
      viscos= VISc
      Unorm = BulkVel
    end if

    Re_tau = u_tau / VISc
    Re_h   = Hnorm * Unorm / VISc

!-----Writing cutline values-------------------------------------------!

!--write file header
!    if(dir == Xdir) XYZ = 'x'
!    if(dir == Ydir) XYZ = 'y'
!    if(dir == Zdir) XYZ = 'z'
!    write(100+Cln,'(A10,I3,A2,A1,A14,F7.3,F7.3,A10,F7.3,F7.3,A2)')       &
!    ' # Cutline', Cln,'; ', XYZ, '-axis (const: ', Jcut, Kcut,           &
!    ' | bounds: ', Icut(1), Icut(1),')'
!
!    if(CLnorm /= PlusM) then
!
!      if(CLnorm==BulkM) write(100+Cln,'(A,1F10.5,A,1F10.5,A)')           &
!      ' # Non-normalized (BulkVel= ', BulkVel, '; Hnorm= ', Hnorm,')'
!
!      if(CLnorm==UtauM) write(100+Cln,'(A,1F10.3,A,1F11.7,A,1PE14.7,A)') &
!      ' # U_tau normalized (Re_tau= ', Re_tau, '; Unorm= ', Unorm,      &
!      '; visc= ', VISc,')'
!
!      if(CLnorm==UtauM .and. u_tau==0.0) write(100+Cln,*)                &
!      '# Normalization value set to 1.0 (u_tau=0.0)!!'
!
!      if(CLnorm==BulkM) write(100+Cln,'(A,1F10.3,A,1F11.7,A,1PE14.7,A)') &
!      ' # Bulk normalized (Re_h= ', Re_h, '; Unorm= ', Unorm,           &
!      '; visc= ', VISc,')'
!
!      write (100+Cln,*) '#'
!      write (100+Cln,'(A3,A1,A39)') ' # ', XYZ,                          &
!      ', U, V, W, T, Kin, Eps, Pk, uv, v_2, f22'
!    else
!      write(100+Cln,'(A,1F10.3,A,1F11.7,A,1PE14.7,A)')                   &
!      ' # Specified normalization (Re_h= ', Re_h,                       &
!      '; Re_tau= ', Re_tau,')'
!
!      write (100+Cln,*) '#'
!      write (100+Cln,'(A3,A1,A23)') ' # ', XYZ, ', Cf, Y+, U/Ub, V/Ub, St'
!    end if

!--browse through cells
    do cp=1,NC                  ! cp - Cutline Point in current CL
      if( Here(cp) ) then

!--parameters
        Tinf  = 20.0
        Flog  = Elog * exp( 9.24 * ((Pr/Prt)**0.75 - 1.0)           &
                        * (1.0 + 0.28*exp(-0.007*Pr/Prt)) * kappa )

!--standard:
        if(CLnorm==UtauM .or. CLnorm==RealM .or. CLnorm==BulkM) then

!--values
!         Yplus = Cmu25 * sqrt(Kin_cl(1)) * WallDs(Ccut(1)) / VISc
!         Twall = Tcl(1) + ( ( Prt / kappa * log(Flog*Yplus) ) * Qflux                       &
!                            / (Cmu25*sqrt(Kin_cl(1))*DENc(material(Nwc))*CAPc(material(Nwc))) )
!         Twall = Tcl(1) + Prt / CAPc(material(1)) * Qflux * WallDs(1) / CONwall(1)

          write(100+Cln, '(2F13.7, 2F7.2, 8F13.7)') Icut(cp)/Hnorm,     &
          Ucl(cp)/Unorm,     Vcl(cp)/Unorm,      VISt_cl(cp)/VISc,     & !Wcl(cp)/Unorm,        &
          ( Tinf - Tcl(cp) ) / Tnorm,                                  &
          Kin_cl(cp)/Unorm**2.0,        Eps_cl(cp)*viscos/Unorm**4.0,  &
          Pk_cl(cp)*viscos/Unorm**4.0,  UiUj_cl(cp)/Unorm**2.0,        &
          v_2_cl(cp)/Unorm**2.0,        f22_cl(cp)*viscos/Unorm**2.0

!--additional:
        else if(CLnorm==PlusM) then

!--values
          Yplus = Cmu25 * sqrt(Kin_cl(cp)) * WallDs(Nwc) / VISc
!          Yplus = ( abs(Ucl(cp)) * WallDs(Nwc) / VISc )**0.5
          Twall = Tcl(cp) + Prt / CAPc(material(Nwc)) * Qflux * WallDs(Nwc) / CONwall(Nwc)
!          Twall = Tcl(cp) + ( ( Prt / kappa * log(Flog*Yplus) ) * Qflux                       &
!                             / (Cmu25*sqrt(Kin_cl(cp))*DENc(material(Nwc))*CAPc(material(Nwc))) )
          if(Yplus > 11.0) then
!            if(WFType==ANL) then
!              TauWall = VISwall(Nwc) * Utan(Nwc) / WallDs(Nwc)
!            else
              TauWall = DENc(material(Nwc)) * kappa           &
                      * sqrt(Kin%n(Nwc)) * Cmu25 * Utan(Nwc)  &
                      / (log(Elog*Yplus))
!            end if
          else
            TauWall = VISc * Utan(Nwc) / WallDs(Nwc)
          end if

          write(100+Cln, '(9F15.7)') Icut(cp)/Hnorm,                         &
          2.0*TauWall/Unorm**2.0, Yplus,                                    & !BulkVel
!          Ucl(cp)/Unorm, Vcl(cp)/Unorm,                                    &
          Kin_cl(cp)/Unorm**2.0, Eps_cl(cp)*viscos/Unorm**4.0,              &
          Pk_cl(cp)*viscos/Unorm**4.0,                                      &
          Qflux/(Twall-Tinf)/(Unorm*DENc(material(Nwc))*CAPc(material(Nwc))), &
          Hnorm / 2.2E-5 * Qflux / ( Twall - Tinf )

        end if              !end if CLmode
      end if               !end if(Here(c))
    end do                !end do cp=1,NC   cutline points

!----------------------------------------------------------------------!

  close(100+Cln)          !close cutline files

  end do                  !end Cln=1,Ncl  number of cutlines
!======================================================================!



!--Go back in T-Rex.cmm file-------------------------------------------!
  if(PRESENT(namAut)) then
    do j=1,Ncl+50
      backspace(7)
    end do
5   CONTINUE
    call ReadC(7,inp,tn,ts,te)
    if(inp(1:len_trim(inp))==ReadLine1(1:len_trim(inp))) then
      backspace(7)
    else
      GOTO 5
    end if
  end if
!----------------------------------------------------------------------!


  RETURN
  END SUBROUTINE CutLine
