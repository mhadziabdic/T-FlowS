!======================================================================!
  subroutine UserCutLines_annulus(grid, name_out) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous direction directions.            
! The results are writen in files name_res.dat and name_res_plus.dat
!----------------------------------------------------------------------!
    use all_mod
    use allp_mod
    use les_mod
    use pro_mod
    use par_mod
    use rans_mod
    use Grid_Mod
!----------------------------------------------------------------------!
    implicit none
!-----------------------------[Arguments]------------------------------!
    type(Grid_Type) :: grid
    real :: Ufric
!------------------------------[Calling]-------------------------------!
    interface
      logical function Approx(A,B,tol)
        real           :: A,B
        real, optional :: tol
      end function Approx
    end interface 
!-------------------------------[Locals]-------------------------------!
    integer             :: Nprob, pl, c, i, count, kk
    character           :: namCoo*80, name_res*80
    character           :: name_res_plus*80
    character, optional :: name_out*(*)
    real,allocatable    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),       & 
                                 var_1(:), var_2(:),     &
                                 var_3(:), Wall_p(:), var_4(:),&
                                 Ufric_p(:), var_5(:), var_6(:),&
                                 var_7(:), var_8(:), var_9(:)
    integer,allocatable :: Np(:), Ncount(:)
    real                :: R, Lscale, R_max
    real                :: Vwall
!======================================================================!

  if(PRESENT(name_out)) then
!---- save the name
    name_res = name_out
    name_res(len_trim(name_out)+1:len_trim(name_out)+8) = '_res.dat'
    name_res_plus = name_out 
    name_res_plus(len_trim(name_out)+1:len_trim(name_out)+13) = '_res_plus.dat'
  else
    name_res = name
    name_res(len_trim(name)+1:len_trim(name)+8) = '_res.dat'
    name_res_plus = name
    name_res_plus(len_trim(name)+1:len_trim(name)+13) = '_res_plus.dat'
  end if
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this_proc < 2) write(6, *) '# Now reading the file:', namCoo
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

    call SSORT (z_p, ind, Nprob, 0)
    
    call CalcShear(U % mean, V % mean, W % mean, Shear)


    allocate(Np(Nprob));    Np=0 
    allocate(Wall_p(Nprob));Wall_p=0.0
    allocate(Ump(Nprob));   Ump=0.0
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(uwp(Nprob));   uwp=0.0
    allocate(vwp(Nprob));   vwp=0.0
    allocate(Ksgsp(Nprob)); Ksgsp=0.0
    allocate(var_1(Nprob));  var_1=0.0
    allocate(var_2(Nprob));  var_2=0.0
    allocate(var_3(Nprob));  var_3=0.0
    allocate(var_4(Nprob));  var_4=0.0
    allocate(var_5(Nprob));  var_5=0.0
    allocate(var_6(Nprob));  var_6=0.0
    allocate(var_7(Nprob));  var_7=0.0
    allocate(var_8(Nprob));  var_8=0.0
    allocate(var_9(Nprob));  var_9=0.0
    allocate(Ufric_p(Nprob));  Ufric_p=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  
    kk = 0

    Lscale = 0.0 
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!

!    do c=-1,-NbC,-1 
!      if(bcmark(c) == 1) then
!        R = sqrt(grid % xc(c)**2+grid % yc(c)**2)
!        Vwall = (-U % n(c) * grid % yc(c) / R  + V % n(c) * grid % xc(c) / R)
!      end if 
!    end do

    Vwall = 1.0

    do c =  1, NC
      R  = (grid % xc(c)*grid % xc(c) + grid % yc(c)*grid % yc(c))**0.5 + tiny
      PHIx(c)   = W % mean(c)
      PHIy(c)   = (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R)
    end do
    call GraPhi(PHIx, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(PHIx, 2, Uy,.TRUE.)    ! dU/dx
    call GraPhi(PHIy, 1, Vx,.TRUE.)    ! dU/dx
    call GraPhi(PHIy, 2, Vy,.TRUE.)    ! dU/dx
    do c =  1, NC
      R  = (grid % xc(c)*grid % xc(c) + grid % yc(c)*grid % yc(c))**0.5 + tiny
      PHIx(c)   = (Ux(c) * grid % xc(c) / R  + Uy(c) * grid % yc(c) / R)
      PHIy(c)   = (Vx(c) * grid % xc(c) / R  + Vy(c) * grid % yc(c) / R)
    end do

    do i = 1, Nprob-1
      do c=1, NC
        Lscale = max(Lscale,WallDs(c))
        R = sqrt(grid % xc(c)**2+grid % yc(c)**2)
           print *, R, z_p(i+1), z_p(i)
        if(R > abs(z_p(i+1)) .and. R < abs(z_p(i))) then
          Wall_p(i) = Wall_p(i) + R 
          if(SIMULA==LES) then
            Ump(i)   = Ump(i) + U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R
            Vmp(i)   = Vmp(i) + (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R)
            Wmp(i)   = Wmp(i) + W % mean(c)
            var_1(i)   = var_1(i) + U % n(c) * grid % xc(c) / R  + V % n(c) * grid % yc(c) / R
            var_2(i)   = var_2(i) + (-U % n(c) * grid % yc(c) / R  + V % n(c) * grid % xc(c) / R)
            var_3(i)   = var_3(i) + W % n(c)
            var_4(i)   = var_4(i) + VISt(c)/VISc
            var_5(i)   = var_5(i) + Ksgs(c)/VISc
            var_6(i)   = var_6(i) + atan(PHIy(c)/PHIx(c))
            uup(i)   = uup(i) + (uu % mean(c)- (U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R)**2.0)
            vvp(i)   = vvp(i) + (vv % mean(c)- (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R)**2.0)
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- (U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R)* &
                                (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R)   )
            uwp(i)   = uwp(i) + (uw % mean(c)- (U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R) * W % mean(c))
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + (VISc * (U % mean(c)**2 + V % mean(c)**2 + W % mean(c)**2)**0.5/WallDs(c))**0.5
            end if
          end if 
 
          if(HOT==YES) then
            if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - (U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            else
              Tmp(i)   = Tmp(i) + T % n(c)
            end if
          end if
          Ncount(i) = Ncount(i) + 1
        end if
      end do 
    end do 

!---- average over all processors
    do pl=1, Nprob-1
      call IGlSum(Ncount(pl))

      call GloSum(Wall_p(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
      call GloSum(var_1(pl))
      call GloSum(var_2(pl))
      call GloSum(var_3(pl))
      call GloSum(var_4(pl))
      call GloSum(var_5(pl))
      call GloSum(var_6(pl))
      call GloSum(Ufric_p(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

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
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        var_4(i)  =  var_4(i)/Ncount(i)
        var_5(i)  =  var_5(i)/Ncount(i)
        var_6(i)  =  var_6(i)/Ncount(i)
        Ufric_p(i)  =  Ufric_p(i)/Ncount(i)
        if(HOT == YES) then
          Tmp(i)    =  Tmp(i)/Ncount(i)
          TTp(i)    =  TTp(i)/Ncount(i)
          uTp(i)    =  uTp(i)/Ncount(i)
          vTp(i)    =  vTp(i)/Ncount(i)
          wTp(i)    =  wTp(i)/Ncount(i)
        end if
      end if
    end do 

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)       
    end do

!    open(12, file='inflow.dat')
!      do i = 1, Nprob-1
!        if(Ncount(i) /= 0) then
!          write(12,'(4E15.7)') var_4(i), var_1(i), var_2(i), var_3(i)
!        end if
!      end do
!    close(12)
 
    open(3,file=name_res)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric*R_max/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(19E15.7)') (Wall_p(i)-0.5)/0.5, Ump(i), Vmp(i)/Vwall, Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), 0.5*(uup(i) + vvp(i) + wwp(i)), var_1(i), var_2(i)/Vwall, var_3(i), var_4(i), var_5(i), var_6(i),&
          atan(2.0*vwp(i)/(wwp(i)-vvp(i))), atan(uvp(i)/uwp(i))
        end if
      end do 
    else if(SIMULA == K_EPS) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:zeta, 8:f' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(12E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i), vwp(i), var_1(i), var_2(i)
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:zeta, 8:f' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    else if(SIMULA == K_EPS_VV) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:v2, 8:f  ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
!          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          write(3,'(9E15.7)') 1.0 - Wall_p(i), abs(Ump(i)), abs(Vmp(i)), -abs(Wmp(i)), 20.0,  &
          uup(i), vvp(i), wwp(i), uvp(i) !, uwp(i)
        end if
      end do 
    else if(SIMULA == K_EPS) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps             ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    else if(SIMULA == SPA_ALL) then
      write(3,'(A1,1X,A50)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:VISt                   ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') Wall_p(i), abs(Ump(i)), abs(Vmp(i)), abs(Wmp(i)),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i)
        end if
      end do 
    end if
    close(3)

    open(3,file=name_res_plus)
    write(3,'(A1,2(A10, F13.5))') '#', 'Utau = ', Ufric, 'Re_tau = ', Ufric*R_max/VISc
    if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA==DNS) then
      write(3,'(A1,2X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:uu+, 6:vv+, 7:ww+, 8:uv+, 9:uw+, 10:vw+, 11:Kin+' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') Wall_p(i)*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)/Ufric**2.0, wwp(i)/Ufric**2.0, uvp(i)/Ufric**2.0, uwp(i)/Ufric**2.0,      &
          vwp(i)/Ufric**2.0, (uup(i) + vvp(i) + wwp(i))/Ufric**2.0
        end if
      end do 
    else if(SIMULA == K_EPS) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:zeta+, 8:f+                     ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(11E15.7)') (Wall_p(i))*Ufric/VISc, Ump(i)/Ufric, Vmp(i)/Ufric, Wmp(i)/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)/Ufric**2.0, wwp(i)/Ufric**2.0, uvp(i)/Ufric**2.0, uwp(i)/Ufric**2.0, &
          vwp(i)/Ufric**2.0, 0.5*(uup(i)/Ufric**2.0+vvp(i)/Ufric**2.0+wwp(i)/Ufric**2.0)
        end if
      end do 
    else if(SIMULA == ZETA) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:zeta+, 8:f+                     ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i), uvp(i)*VISc/Ufric**2.0, uwp(i)/VISc
        end if
      end do 
    else if(SIMULA == K_EPS_VV) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:Kin+, 6:Eps+, 7:v2+, 8:f+                       ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i)/Ufric**2.0, uvp(i)*VISc/Ufric**2.0, uwp(i)/Ufric**2.0
        end if
      end do 
    else if(SIMULA == SPA_ALL) then
      write(3,'(A1,1X,A80)') '#', '1:Xrad+, 2:U+, 3:V+, 4:W+, 5:VISt/VISc                                       ' 
      do i = 1, Nprob-1
        if(Ncount(i) /= 0) then
          write(3,'(9E15.7)') (Wall_p(i))*Ufric/VISc, abs(Ump(i))/Ufric, abs(Vmp(i))/Ufric, abs(Wmp(i))/Ufric,   &
          uup(i)/Ufric**2.0, vvp(i)*VISc/Ufric**4.0, wwp(i), uvp(i)*VISc/Ufric**2.0, uwp(i)/Ufric**2.0
        end if
      end do 
    end if
    close(3)

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
    deallocate(var_1)
    deallocate(var_2)
    deallocate(var_3)
    deallocate(var_4)
    deallocate(var_5)
    deallocate(var_6)
    deallocate(Ksgsp)
    deallocate(Ufric_p)
    deallocate(Wall_p)
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  end subroutine UserCutLines_annulus
