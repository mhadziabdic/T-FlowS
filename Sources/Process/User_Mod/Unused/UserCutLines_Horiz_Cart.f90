!======================================================================!
  subroutine UserCutLines_Horiz_Cart(grid) 
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
!-----------------------------[Arguments]------------------------------!
!  real :: y(-NbC:NC)
  real :: Rad_2, Ufric 
!------------------------------[Calling]-------------------------------!
  integer             :: Nprob, pl, c, i, count, k, c1, c2, s, N_hor
  character           :: namCoo*80, namPro*80, answer*80, JetIn*14, namOut*80
  real,allocatable    :: r1_p(:), r2_p(:), z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_mp(:), &
                                 var_4(:), var_5(:)  
  integer,allocatable :: Np(:), Ncount(:)
  real                :: R, Urad_n, Utan_n, R1, R2, Urad, Utan, dummy 
!======================================================================!

  open(9, file='Horiz_positions.dat')
  if(this < 2) print *, '# Now reading the file: Horiz_positions.dat ' 
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
    if(this < 2) print *, '# Now reading the file: rad_coordinate.dat ' 
    open(9, file='Z_coordinate.dat')

!---- write the number of probes 
    read(9,*) Nprob
    allocate(z_p(Nprob))

!---- write the probe coordinates out
    do pl=1,Nprob
      read(9,*) dummy, z_p(pl)
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
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do k = 1, N_hor
    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = zc(c) 
        if(xc(c) < r1_p(k) .and. xc(c) > r2_p(k)) then
          if(Rad_2 > z_p(i) .and. Rad_2 < z_p(i+1)) then
 
            Ump(i)   = Ump(i) + U%n(c) 
            Vmp(i)   = Vmp(i) + V%n(c)
            Wmp(i)   = Wmp(i) + W%n(c)
            if(SIMULA==EBM.or.SIMULA==HJ)then
              uuP(i)   = uup(i) + uu%n(c) 
              vvp(i)   = vvp(i) + vv%n(c)
              wwp(i)   = wwp(i) + ww%n(c)
              uvp(i)   = uvp(i) + uv%n(c)
              uwp(i)   = uwp(i) + uw%n(c)
              var_1(i) = var_1(i) + Kin%n(c)
              var_2(i) = var_2(i) + Eps%n(c) 
              var_3(i) = var_3(i) + VISt(c)/VISc
            else if(SIMULA==ZETA) then 
              uup(i)   = uup(i) + 2.0*VISt(c)*Ux(c) 
              vvp(i)   = vvp(i) + 2.0*VISt(c)*Vy(c)
              wwp(i)   = wwp(i) + 2.0*VISt(c)*Wz(c)
              uvp(i)   = uvp(i) + VISt(c)*(Uy(c) + Vx(c))
              uwp(i)   = uwp(i) + VISt(c)*(Uz(c) + Wx(c))
              vwp(i)   = vwp(i) + VISt(c)*(Vz(c) + Wy(c))
              var_1(i) = var_1(i) + Kin%n(c)
              var_2(i) = var_2(i) + Eps%n(c) 
              var_3(i) = var_3(i) + VISt(c)/VISc
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
      JetIn  = 'BSt_Hor_x'
      write(JetIn(9:13),'(I1,A4)') k, '.dat'
    else
      JetIn  = 'BSt_Hor_xx'
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
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        var_4(i)  =  var_4(i)/Ncount(i)
        var_5(i)  =  var_5(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)

        write(3,'(12E15.7, I7)') (z_p(i)+z_p(i+1))/2.0, Ump(i), Vmp(i), Wmp(i) , uup(i), vvp(i), &
                                 wwp(i), uvp(i), uwp(i), var_1(i), var_2(i), var_3(i), Ncount(i)

        Wmp(i)    = 0.0   
        Ump(i)    = 0.0 
        Vmp(i)    = 0.0 
        uup(i)    = 0.0 
        vvp(i)    = 0.0 
        wwp(i)    = 0.0 
        uvp(i)    = 0.0 
        uwp(i)    = 0.0 
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


  if(this < 2) print *, 'Finished with UserCutLines_Horiz_Cart'

  end subroutine UserCutLines_Horiz_Cart
