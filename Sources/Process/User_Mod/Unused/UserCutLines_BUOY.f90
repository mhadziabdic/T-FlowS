!======================================================================!
  subroutine UserCutLines_BUOY(grid, y) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous directions. 
! The results are writen in files name_buoy_res.dat and name_buoy_res_plus.dat
!----------------------------------------------------------------------!
    use all_mod
    use allp_mod
    use les_mod
    use pro_mod
    use par_mod
    use rans_mod
!----------------------------------------------------------------------!
    implicit none
!-----------------------------[Arguments]------------------------------!
    real :: y(-NbC:NC)
    real :: Ufric, Wall_near 
!------------------------------[Calling]-------------------------------!
    interface
      logical function Approx(A,B,tol)
        real           :: A,B
        real, optional :: tol
      end function Approx
    end interface 
!-------------------------------[Locals]-------------------------------!
    integer             :: Nprob, pl, c, i, count, s, c1, c2
    character           :: namCoo*80, namPro*80, answer*80, namRes*80
    character           :: namRes_plus*80
    real,allocatable    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 uupR(:), vvpR(:), wwpR(:), &
                                 uvpR(:), uwpR(:), vwpR(:), &
                                 Tmp(:), TTmp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 uTpR(:), vTpR(:), wTpR(:), &
                                 Ksgsp(:), ind(:),       & 
                                 VIStp(:), var_2(:), var_3(:), &
                                 Wall_p(:), Ufric_p(:), Kinp(:), Epsp(:), &
                                 v2p(:), fp(:) 
    integer,allocatable :: Np(:), Ncount(:)
    real                :: R, Urad_mean, Utan_mean, dummy, Lscale, Twall
!======================================================================!

    namPro = name

    namRes = name
    namRes(len_trim(name)+1:len_trim(name)+13) = '_buoy_res.dat'
    namRes_plus = name
    namRes_plus(len_trim(name)+1:len_trim(name)+18) = '_buoy_res_plus.dat'
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this < 2)  print *, '# Now reading the file:', namCoo
    open(9, file=namCoo)

!---- write the number of searching intervals 
    read(9,*) Nprob
    allocate(z_p(Nprob*2))
    allocate(ind(Nprob*2))

!---- read the intervals positions
    do pl=1,Nprob
      read(9,*) ind(pl), z_p(pl) 
    end do
    close(9)

    call GraPhi(U % n, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(U % n, 2, Uy,.TRUE.)    ! dU/dy
    call GraPhi(U % n, 3, Uz,.TRUE.)    ! dU/dz
    call GraPhi(V % n, 1, Vx,.TRUE.)    ! dV/dx
    call GraPhi(V % n, 2, Vy,.TRUE.)    ! dV/dy
    call GraPhi(V % n, 3, Vz,.TRUE.)    ! dV/dz
    call GraPhi(W % n, 1, Wx,.TRUE.)    ! dW/dx
    call GraPhi(W % n, 2, Wy,.TRUE.)    ! dW/dy
    call GraPhi(W % n, 3, Wz,.TRUE.)    ! dW/dz

    call GraPhi(T % n,1,PHIx,.TRUE.)
    call GraPhi(T % n,2,PHIy,.TRUE.)
    call GraPhi(T % n,3,PHIz,.TRUE.)

    call SSORT (z_p, ind, Nprob, 0)
    print *, Nprob
    
    call CalcShear(U % n, V % n, W % n, Shear)

    if(SIMULA==EBM.and.MODE/=HYB) then
      do c=1,NC
        VAR1x(c) = 0.22*Tsc(c) *&
                   (uu%n(c)*PHIx(c)+uv%n(c)*PHIy(c)+uw%n(c)*PHIz(c))

        VAR1y(c) = 0.22*Tsc(c)*&
                   (uv%n(c)*PHIx(c)+vv%n(c)*PHIy(c)+vw%n(c)*PHIz(c))

        VAR1z(c) = 0.22*Tsc(c)*&
                   (uw%n(c)*PHIx(c)+vw%n(c)*PHIy(c)+ww%n(c)*PHIz(c))
      end do
 
      call GraPhi(VAR1x,1,VAR2x,.TRUE.)
      call GraPhi(VAR1y,2,VAR2y,.TRUE.)
      call GraPhi(VAR1z,3,VAR2z,.TRUE.)
    end if

    allocate(Np(Nprob));      Np=0 
    allocate(Wall_p(Nprob));  Wall_p=0.0
    allocate(Ump(Nprob));     Ump=0.0
    allocate(Vmp(Nprob));     Vmp=0.0
    allocate(Wmp(Nprob));     Wmp=0.0
    allocate(uup(Nprob));     uup=0.0
    allocate(vvp(Nprob));     vvp=0.0
    allocate(wwp(Nprob));     wwp=0.0
    allocate(uvp(Nprob));     uvp=0.0
    allocate(uwp(Nprob));     uwp=0.0
    allocate(vwp(Nprob));     vwp=0.0
    allocate(uupR(Nprob)); uup=0.0
    allocate(vvpR(Nprob)); vvp=0.0
    allocate(wwpR(Nprob)); wwp=0.0
    allocate(uvpR(Nprob)); uvp=0.0
    allocate(uwpR(Nprob)); uwp=0.0
    allocate(vwpR(Nprob)); vwp=0.0
    allocate(Kinp(Nprob));    Kinp=0.0
    allocate(Epsp(Nprob));    Epsp=0.0
    allocate(v2p(Nprob));     v2p=0.0
    allocate(fp(Nprob));      fp=0.0
    allocate(Ksgsp(Nprob));   Ksgsp=0.0
    allocate(VIStp(Nprob));   VIStp=0.0
    allocate(Ufric_p(Nprob)); Ufric_p = 0.0
    allocate(Ncount(Nprob));  Ncount=0
    allocate(Tmp(Nprob));     Tmp=0.0
    allocate(TTmp(Nprob));     TTmp=0.0
    allocate(uTp(Nprob));     uTp=0.0
    allocate(vTp(Nprob));     vTp=0.0
    allocate(wTp(Nprob));     wTp=0.0
    allocate(uTpR(Nprob)); uTpR=0.0
    allocate(vTpR(Nprob)); vTpR=0.0
    allocate(wTpR(Nprob)); wTpR=0.0
    
    count = 0
    Lscale = 0.0
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
    do i = 1, Nprob-1
      do c=1, NC
        Lscale = max(Lscale,WallDs(c))
        if(y(c) > (z_p(i+1)) .and. y(c) < (z_p(i))) then
          Wall_p(i) = Wall_p(i) + WallDs(c) 
          Ump(i)   = Ump(i) + U % n(c)
          Vmp(i)   = Vmp(i) + V % n(c)
          Wmp(i)   = Wmp(i) + W % n(c)
          Kinp(i)  = Kinp(i) + Kin % n(c) 
          Epsp(i)  = Epsp(i) + Eps % n(c) 
          v2p(i)   = v2p(i)  + v_2 %n(c) !- CONc(material(c))*PHIz(c)  !v_2 % n(c) 
          fp(i)    = fp(i)  + f22%n(c)   !tt % n(c) ! - VISt(c)*PHIz(c)
!          v2p(i)   = v2p(i) + v_2 % n(c) 
!          fp(i)    = fp(i) + VISt(c)/VISc !f22 % n(c)
          uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
          vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
          wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
          uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
          uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
          vwp(i)   = vwp(i) + VISt(c)*(Vz(c) + Wy(c))
          uupR(i)  = uupR(i) + uu_res(c) 
          vvpR(i)  = vvpR(i) + vv_res(c)
          wwpR(i)  = wwpR(i) + ww_res(c)
          uvpR(i)  = uvpR(i) + uv_res(c)
          uwpR(i)  = uwpR(i) + uw_res(c)
          vwpR(i)  = vwpR(i) + vw_res(c)
          VIStp(i) = VIStp(i) + VISt(c)/VISc 
          if(IsNearWall(c)) then
            if(Ynd(c)>5.0) then      
              Ufric_p(i) = Ufric_p(i) + sqrt(TauWall(c))
            else
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % n(c)**2 + V % n(c)**2 + &
              W% n(c)**2)**0.5/WallDs(c))**0.5      
            end if
          end if 
          Tmp(i)   = Tmp(i) + T % n(c)
          TTmp(i)   = TTmp(i) + TT % n(c)
          uTp(i)   = uTp(i) + uT % n(c)
          vTp(i)   = vTp(i) + vT % n(c)
          wTp(i)   = wTp(i) + wT % n(c)
          uTpR(i)   = uTpR(i) + uT_res(c)
          vTpR(i)   = vTpR(i) + vT_res(c)
          wTpR(i)   = wTpR(i) + wT_res(c)
          Ncount(i) = Ncount(i) + 1
        end if
      end do 
    end do 
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2  < 0) then
        if( TypeBC(c2) ==  WALL.or. TypeBC(c2) ==  WALLFL) then
          Twall = T % n(c2) 
        end if
      end if
    end do

!---- average over all processors
    do pl=1, Nprob-1
      call IGlSum(Ncount(pl))

      call GloSum(Wall_p(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))

      call GloSum(Kinp(pl))
      call GloSum(Epsp(pl))
      call GloSum(fp(pl))
      call GloSum(v2p(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))
      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))

      call GloSum(uupR(pl))
      call GloSum(vvpR(pl))
      call GloSum(wwpR(pl))
      call GloSum(uvpR(pl))
      call GloSum(uwpR(pl))
      call GloSum(vwpR(pl))

      call GloSum(VIStp(pl))
      call GloSum(Ufric_p(pl))

      call GloSum(Tmp(pl))
      call GloSum(TTmp(pl))
      call GloSum(uTp(pl))
      call GloSum(vTp(pl))
      call GloSum(wTp(pl))
      call GloSum(uTpR(pl))
      call GloSum(vTpR(pl))
      call GloSum(wTpR(pl))

      count =  count + Ncount(pl) 
    end do

    call wait

    do i = 1, Nprob-1
      if(Ncount(i) /= 0) then
        Wall_p(i) =  Wall_p(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        Wmp(i)    =  Wmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        vwp(i)    =  vwp(i)/Ncount(i)
        uupR(i)   =  uupR(i)/Ncount(i)
        vvpR(i)   =  vvpR(i)/Ncount(i)
        wwpR(i)   =  wwpR(i)/Ncount(i)
        uvpR(i)   =  uvpR(i)/Ncount(i)
        uwpR(i)   =  uwpR(i)/Ncount(i)
        vwpR(i)   =  vwpR(i)/Ncount(i)
        Kinp(i)   =  Kinp(i)/Ncount(i)
        Epsp(i)   =  Epsp(i)/Ncount(i)
        v2p(i)    =  v2p(i)/Ncount(i)
        fp(i)     =  fp(i)/Ncount(i)
        VIStp(i)  =  VIStp(i)/Ncount(i)
        Ufric_p(i)=  Ufric_p(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        TTmp(i)   =  TTmp(i)/Ncount(i)
        uTp(i)    =  uTp(i)/Ncount(i)
        vTp(i)    =  vTp(i)/Ncount(i)
        wTp(i)    =  wTp(i)/Ncount(i)
        uTpR(i)   =  uTpR(i)/Ncount(i)
        vTpR(i)   =  vTpR(i)/Ncount(i)
        wTpR(i)   =  wTpR(i)/Ncount(i)
      end if
    end do 

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)
    end do
10  continue

    if(Ufric == 0.0) then
      if(this < 2) print *, 'Friction velocity is zero in UserCutLines_BUOY.f90 !'
      return
    end if

    open(3,file=namRes)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    write(3,'(A1,1X,A140)') '#', 'Xrad U V W T (-5) TT uu vv ww uv (-10) &
                             uw vw uuR vvR wwR (-15) uvR uwR vwR Kin Eps (-20) &
                             v2 f22 VISratio uT vT (-25) wT uTR vTR wTR (-29)'
    do i = Nprob, 1, -1 
      if(Ncount(i) /= 0) then
        write(3,'(29E16.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), TTmp(i), &
        uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), uupR(i), vvpR(i), wwpR(i), &
        uvpR(i), uwpR(i), vwpR(i), Kinp(i), Epsp(i), v2p(i), fp(i), VIStp(i), uTp(i), &
        vTp(i), wTp(i), uTpR(i), vTpR(i), wTpR(i)
      end if
    end do     
    close(3)

    open(3,file=namRes_plus)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric/VISc
    write(3,'(A1,2(A10, F14.8))') '#', 'Tflux = ', Tflux, ' K+ =', Tflux/Ufric * (VISc/CONc(material(1)))**ONE_THIRD
    write(3,'(A1,1X,A130)') '#', 'Xrad U V W T (-5) uu vv ww uv uw (-10) &
                             vw Kin Eps v2 f22 (-15) VISratio uT vT wT uTR (-20) &
                             vTR wTR Tnorm sqr(TT)norm (-25)' 
    do i = Nprob, 1, -1 
      if(Ncount(i) /= 0) then
        write(3,'(25E16.7)') Wall_p(i)*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric, (Twall - Tmp(i))*Ufric/Tflux, &
        uup(i)/Ufric**2, vvp(i)/Ufric**2, wwp(i)/Ufric**2, uvp(i)/Ufric**2, uwp(i)/Ufric**2, vwp(i)/Ufric**2, &
        uupR(i)/Ufric**2, vvpR(i)/Ufric**2, wwpR(i)/Ufric**2, uvpR(i)/Ufric**2, uwpR(i)/Ufric**2, vwpR(i)/Ufric**2, &
        Kinp(i)/Ufric**2, Epsp(i)*VISc/Ufric**4.0, v2p(i), Ufric*fp(i)/Tflux, VIStp(i), &
        uTp(i)/Tflux, vTp(i)/Tflux, wTp(i)/Tflux, uTpR(i)/Tflux, vTpR(i)/Tflux, wTpR(i)/Tflux, &
        (Tmp(i)-95.0)/10.0, sqrt(TTmp(i))/5.0
      end if
    end do 
    
    close(3)

    deallocate(Ufric_p);
    deallocate(Ncount);
    deallocate(Wall_p);
    deallocate(ind)
    deallocate(Np)
    deallocate(z_p)
    deallocate(Ump)
    deallocate(Vmp)
    deallocate(Wmp)
    deallocate(uup)
    deallocate(vvp)
    deallocate(wwp)
    deallocate(uvp)
    deallocate(uwp)
    deallocate(vwp)
    deallocate(uupR)
    deallocate(vvpR)
    deallocate(wwpR)
    deallocate(uvpR)
    deallocate(uwpR)
    deallocate(vwpR)
    deallocate(Kinp)
    deallocate(Epsp)
    deallocate(v2p)
    deallocate(fp)
    deallocate(VIStp)
    deallocate(Ksgsp)
    deallocate(Tmp)
    deallocate(TTmp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
    deallocate(uTpR)
    deallocate(vTpR)
    deallocate(wTpR)
    
  end subroutine UserCutLines_BUOY
