!======================================================================!
  subroutine UserCutLines_Y_dir() 
!----------------------------------------------------------------------!
! Reads the ".1D" file created by the "Generator" and averages the     !
! results in the planes defined by coordinates in it. Then averages    !
! the values of Un, Vn, Wn, uu, vv, ww, uv, uw and vw and     !
! writes them into file ".1Dr".                                        !
!----------------------------------------------------------------------!
  use all_mod
  use allp_mod
  use les_mod
  use pro_mod
  use par_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
!  real :: y(-NbC:NC)
  real :: Rad_2, Ufric 
!------------------------------[Calling]-------------------------------!
  integer             :: Nprob, pl, c, dummy, i, count, k, c1, c2, s, N_hor
  character           :: namCoo*80, namPro*80, answer*80, JetIn*14, namOut*80
  real,allocatable    :: r1_p(:), r2_p(:), z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 uum(:), vvm(:), wwm(:), &
                                 uvm(:), uwm(:), vwm(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_mp(:), &
                                 var_4(:), var_5(:)  
  integer,allocatable :: Np(:), Ncount(:)
  real                :: R, Urad_n, Utan_n, R1, R2, Urad, Utan 
!======================================================================!

  open(9, file='Hor_positions.dat')
  if(this < 2) write(6, *) '# Now reading the file: Hor_positions.dat ' 
  read(9,*) N_hor
  allocate(r1_p(N_hor))
  allocate(r2_p(N_hor))
  do pl=1, N_hor
    read(9,*) r1_p(pl), r2_p(pl)
  end do
  close(9)

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    if(this < 2) write(6, *) '# Now reading the file: Y_coordinate.dat ' 
    open(9, file='Y_coordinate.dat')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) z_p(pl)
    end do
    close(9)

    allocate(Np(Nprob));    Np=0 
    allocate(Ump(Nprob));   Ump=0.0
    allocate(Vmp(Nprob));   Vmp=0.0
    allocate(Wmp(Nprob));   Wmp=0.0
    allocate(uup(Nprob));   uup=0.0
    allocate(vvp(Nprob));   vvp=0.0
    allocate(wwp(Nprob));   wwp=0.0
    allocate(uvp(Nprob));   uvp=0.0
    allocate(uwp(Nprob));   uwp=0.0
    allocate(vwp(Nprob));   vwp=0.0
    allocate(uum(Nprob));   uum=0.0
    allocate(vvm(Nprob));   vvm=0.0
    allocate(wwm(Nprob));   wwm=0.0
    allocate(uvm(Nprob));   uvm=0.0
    allocate(uwm(Nprob));   uwm=0.0
    allocate(vwm(Nprob));   vwm=0.0
    allocate(Ksgsp(Nprob)); Ksgsp=0.0
    allocate(Rad_mp(Nprob));  Rad_mp=0.0
    allocate(var_1(Nprob));  var_1=0.0
    allocate(var_2(Nprob));  var_2=0.0
    allocate(var_3(Nprob));  var_3=0.0
    allocate(var_4(Nprob));  var_4=0.0
    allocate(var_5(Nprob));  var_5=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  
    call CalcShear(U % n, V % n, W % n, Shear)

!=================================================================!
! Part related to V-budget calculation
!=================================================================!
    call GradP(P % mean,Px,Py,Pz)  !IZMENA
    call GraPhi(V % mean, 2, PHI1z,.TRUE.)
    call GraPhi(V % mean, 1, PHI2z,.TRUE.)
    call GraPhi(U % mean, 2, PHI3x,.TRUE.)
    if(SIMULA == LES) then
      do c = 1,NC
        PHI1x(c) = (vv % mean(c)- V % mean(c) * V % mean(c))
        PHI2x(c) = (uv % mean(c)- U % mean(c) * V % mean(c))
      end do
    else if(SIMULA==HJ) then
      do c = 1,NC
        PHI1x(c) = VAR10y(c) - V%mean(c)*V%mean(c) + vv%mean(c)
        PHI2x(c) = VAR11x(c) - U % mean(c) * V % mean(c) + uv%mean(c) 
      end do
    else if(SIMULA==ZETA) then
      do c = 1,NC
        PHI1x(c) = (vv % mean(c)- V % mean(c) * V % mean(c)) + &
        0.66666*Kin%mean(c) - &
        2*CmuD*Kin%mean(c)**2*v_2%mean(c)/Eps%mean(c)*PHI1z(c) 
        PHI2x(c) = (uv % mean(c)- U % mean(c) * V % mean(c)) - &
        CmuD*Kin%mean(c)**2*v_2%mean(c)/Eps%mean(c)*(PHI2z(c)+PHI3x(c))
      end do
    end if        
           
    call GraPhi(PHI1x, 2, PHI1y,.TRUE.)
    call GraPhi(PHI2x, 1, PHI3x,.TRUE.)
!=================================================================!
!=================================================================!
! Part related to U-budget calculation
!=================================================================!
!    call GradP(P % mean,Px,Py,Pz)  !IZMENA
!    call GraPhi(U % mean, 1, PHI1z,.TRUE.)
!    call GraPhi(U % mean, 2, PHI2z,.TRUE.)
!    call GraPhi(U % mean, 2, PHI3x,.TRUE.)
!    if(SIMULA == LES) then
!      do c = 1,NC
!        PHI1x(c) = (uu % mean(c)- U % mean(c) * U % mean(c))
!        PHI2x(c) = (uv % mean(c)- U % mean(c) * V % mean(c))
!      end do
!    else if(SIMULA==HJ) then
!      do c = 1,NC
!        PHI1x(c) = VAR10x(c) - U%mean(c)*U%mean(c) + uu%mean(c)
!        PHI2x(c) = VAR11x(c) - U % mean(c) * V % mean(c) + uv%mean(c) 
!      end do
!    else if(SIMULA==ZETA) then
!      do c = 1,NC
!        PHI1x(c) = (vv % mean(c)- V % mean(c) * V % mean(c)) + &
!        0.66666*Kin%mean(c) - &
!        2*CmuD*Kin%mean(c)**2*v_2%mean(c)/Eps%mean(c)*PHI1z(c) 
!        PHI2x(c) = (uv % mean(c)- U % mean(c) * V % mean(c)) - &
!        CmuD*Kin%mean(c)**2*v_2%mean(c)/Eps%mean(c)*(PHI2z(c)+PHI3x(c))
!      end do
!    end if        
!           
!    call GraPhi(PHI1x, 1, PHI1y,.TRUE.)
!    call GraPhi(PHI2x, 2, PHI3x,.TRUE.)
!=================================================================!
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do k = 1, N_hor
    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = yc(c) 
        if(xc(c) < r1_p(k) .and. xc(c) > r2_p(k)) then
          if(Rad_2 > z_p(i) .and. Rad_2 < z_p(i+1)) then
            if(SIMULA==LES) then      
              Ump(i)   = Ump(i) + U % mean(c)
              Vmp(i)   = Vmp(i) + V % mean(c)
              Wmp(i)   = Wmp(i) + W % mean(c)
              uup(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
              vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
              wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
              uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
              uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
              vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))
!              var_1(i) = var_1(i) + V % mean(c)*PHI1z(c) + U%mean(c)*PHI2z(c)
!              var_2(i) = var_2(i) + Py(c)
!              var_3(i) = var_3(i) + PHI1y(c) + PHI3x(c)
              var_1(i) = var_1(i) + U % mean(c)*PHI1z(c) + V%mean(c)*PHI2z(c)
              var_2(i) = var_2(i) + Px(c)
              var_3(i) = var_3(i) + PHI1y(c) + PHI3x(c)
              var_4(i) = var_4(i) + 0.5*((uu%mean(c)-U%mean(c)*U%mean(c)) + &
                                         (vv%mean(c)-V%mean(c)*V%mean(c)) + &
                                         (ww%mean(c)-W%mean(c)*W%mean(c))) 
            else if(SIMULA == HJ) then  
              Ump(i)   = Ump(i) + U % mean(c)
              Vmp(i)   = Vmp(i) + V % mean(c)
              Wmp(i)   = Wmp(i) + W % mean(c)
              uup(i)   = uup(i) + uu % mean(c)
              vvp(i)   = vvp(i) + vv % mean(c)
              wwp(i)   = wwp(i) + ww % mean(c)
              uvp(i)   = uvp(i) + uv % mean(c)
              uwp(i)   = uwp(i) + uw % mean(c)
              vwp(i)   = vwp(i) + vw % mean(c)
              uum(i)   = uum(i) + VAR10x(c) - U % mean(c) * U % mean(c)
              vvm(i)   = vvm(i) + VAR10y(c) - V % mean(c) * V % mean(c)
              wwm(i)   = wwm(i) + VAR10z(c) - W % mean(c) * W % mean(c)
              uvm(i)   = uvm(i) + VAR11x(c) - U % mean(c) * V % mean(c)
              uwm(i)   = uwm(i) + VAR11y(c) - U % mean(c) * W % mean(c)
              vwm(i)   = vwm(i) + VAR11z(c) - V % mean(c) * W % mean(c)

              var_1(i) = var_1(i) + 0.5*(uu%mean(c)+vv%mean(c)+ww%mean(c)) 
              var_2(i) = var_2(i) + V % mean(c)*PHI1z(c) + U%mean(c)*PHI2z(c)
              var_3(i) = var_3(i) + Py(c)
              var_4(i) = var_4(i) + PHI1y(c) + PHI3x(c)

!              var_2(i) = var_2(i) + U%mean(c)*PHI1z(c) + V%mean(c)*PHI2z(c)
!              var_3(i) = var_3(i) + Px(c)
!              var_4(i) = var_4(i) + PHI1y(c) + PHI3x(c)
!              var_3(i) = var_3(i) + VISt(c)/VISc
!              var_3(i) = var_3(i) + Kin % mean(c) + 0.5*((uu%mean(c)-U%mean(c)*U%mean(c)) + &
!                                         (vv%mean(c)-V%mean(c)*V%mean(c)) + &
!                                         (ww%mean(c)-W%mean(c)*W%mean(c))) 
!              var_4(i) = var_4(i) + V % mean(c)*PHI1z(c) + U%mean(c)*PHI2z(c)
!              var_5(i) = var_5(i) + PHI1y(c) + PHI3x(c) 
            end if        
            Rad_mp(i) = Rad_mp(i) + zc(c) 
            Ncount(i) = Ncount(i) + 1
          end if
        end if
      end do 
    end do 

!---- average over all processors
    do pl=1, Nprob
      call IGlSum(Ncount(pl))
      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))
      call GloSum(Rad_mp(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))
      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
      call GloSum(uum(pl))
      call GloSum(vvm(pl))
      call GloSum(wwm(pl))
      call GloSum(uvm(pl))
      call GloSum(uwm(pl))
      call GloSum(vwm(pl))
      call GloSum(var_1(pl))
      call GloSum(var_2(pl))
      call GloSum(var_3(pl))
      call GloSum(var_4(pl))
      call GloSum(var_5(pl))

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do

    if(k < 10) then
      JetIn  = 'Cut_Ver_x'
      write(JetIn(9:13),'(I1,A4)') k, '.dat'
    else
      JetIn  = 'Cut_Ver_xx'
      write(JetIn(9:14),'(I2,A4)') k, '.dat'
    end if 

    open(3,file=JetIn)
    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        vwp(i)    =  vwp(i)/Ncount(i)
        uum(i)    =  uum(i)/Ncount(i)
        vvm(i)    =  vvm(i)/Ncount(i)
        wwm(i)    =  wwm(i)/Ncount(i)
        uvm(i)    =  uvm(i)/Ncount(i)
        uwm(i)    =  uwm(i)/Ncount(i)
        vwm(i)    =  vwm(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        var_4(i)  =  var_4(i)/Ncount(i)
        var_5(i)  =  var_5(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)

        write(3,'(21E15.7, I7)') (z_p(i)+z_p(i+1))/2.0, Ump(i), Vmp(i), Wmp(i), uup(i), vvp(i), &
                                wwp(i), uvp(i), uwp(i), vwp(i), uum(i), &
                                vvm(i),wwm(i), uvm(i), uwm(i), vwm(i), &
                                var_1(i), 0.5*(uum(i)+vvm(i)+wwm(i)), var_2(i), var_3(i), var_4(i), Ncount(i)

        Wmp(i)    = 0.0   
        Ump(i)    = 0.0 
        Vmp(i)    = 0.0 
        uup(i)    = 0.0 
        vvp(i)    = 0.0 
        wwp(i)    = 0.0 
        uvp(i)    = 0.0 
        uwp(i)    = 0.0 
        vwp(i)    = 0.0 
        uum(i)    = 0.0 
        vvm(i)    = 0.0 
        wwm(i)    = 0.0 
        uvm(i)    = 0.0 
        uwm(i)    = 0.0 
        vwm(i)    = 0.0 
        var_1(i)  = 0.0 
        var_2(i)  = 0.0 
        var_3(i)  = 0.0
        var_4(i)  = 0.0
        var_5(i)  = 0.0
        Rad_mp(i) = 0.0 
        Ncount(i) = 0
      end if
    end do 
    close(3)
  end do
 
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
  deallocate(uum)
  deallocate(vvm)
  deallocate(wwm)
  deallocate(uvm)
  deallocate(uwm)
  deallocate(vwm)
  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  deallocate(var_4)
  deallocate(var_5)
  deallocate(Rad_mp)
  deallocate(Ksgsp)
  deallocate(Ncount)
  if(HOT==YES) then
    deallocate(Tmp)
    deallocate(TTp)
    deallocate(uTp)
    deallocate(vTp)
    deallocate(wTp)
  end if


  if(this < 2) write(*,*) 'Finished with UserCutLines_Y_dir'

  end subroutine UserCutLines_Y_dir
