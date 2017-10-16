!======================================================================!
  subroutine UserCutLines_RB_conv() 
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
                                 var_4(:), var_5(:), var_6(:)  
  integer,allocatable :: Np(:), Ncount(:)
  real                :: R, Urad_n, Utan_n, R1, R2, Urad, Utan 
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
    if(this_proc < 2) write(6, *) 'Now reading the file:', namCoo
    open(9, file=namCoo)

!---- write the number of probes 
    read(9,'(I8)') Nprob
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

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  

    do i = 1, Nprob-1
      do c=1,NC
        Rad_2 = zc(c) 
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
            var_1(i) = var_1(i) + 0.5*((uu%mean(c)-U%mean(c)*U%mean(c)) + &
                                       (vv%mean(c)-V%mean(c)*V%mean(c)) + &
                                       (ww%mean(c)-W%mean(c)*W%mean(c))) 
            Tmp(i)   = Tmp(i) + (T % mean(c)-5.0)/10.0
            TTp(i)   = TTp(i) + (tt % mean(c) - T%mean(c)*T % mean(c))
            uTp(i)   = uTp(i) + (ut % mean(c) - U%mean(c)*T % mean(c)) 
            vTp(i)   = vTp(i) + (vt % mean(c) - V%mean(c)*T % mean(c))
            wTp(i)   = wTp(i) + (wt % mean(c) - W%mean(c)*T % mean(c)) 
            Rad_mp(i)= Rad_mp(i) + zc(c) 
            Ncount(i)= Ncount(i) + 1
          else
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
            uup(i)   = uup(i) + Kin %n(c) !(uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + Eps % n(c) !(vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + uw %n(c)   !+ (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + VISt(c)/VISc !+ (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) !+ (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) !+ (vw % mean(c)- V % mean(c) * W % mean(c))
!            var_1(i) = var_1(i) + 0.5*((uu%mean(c)-U%mean(c)*U%mean(c)) + &
!                                       (vv%mean(c)-V%mean(c)*V%mean(c)) + &
!                                       (ww%mean(c)-W%mean(c)*W%mean(c))) 
            Tmp(i)   = Tmp(i) + (T % n(c)-5.0)/10.0
            TTp(i)   = TTp(i) + tt % n(c) !sqrt(tt % n(c))/10.0 !- T%mean(c)*T % mean(c))
            uTp(i)   = uTp(i) + ut % n(c) !- U%mean(c)*T % mean(c)) 
            vTp(i)   = vTp(i) + vt % n(c) !- V%mean(c)*T % mean(c))
            wTp(i)   = wTp(i) + wt % n(c) !- W%mean(c)*T % mean(c)) 
            Rad_mp(i)= Rad_mp(i) + zc(c) 
            Ncount(i)= Ncount(i) + 1
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

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
        call GloSum(TTp(pl))
        call GloSum(uTp(pl))
        call GloSum(vTp(pl))
        call GloSum(wTp(pl))
      end if
    end do
    k = 1
    JetIn  = 'Cut_Ver_Y'
    write(JetIn(9:13),'(I1,A4)') k, '.dat'

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
        var_1(i)  =  var_1(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i) 
        TTp(i)    =  TTp(i)/Ncount(i) 
        uTp(i)    =  uTp(i)/Ncount(i) 
        vTp(i)    =  vTp(i)/Ncount(i) 
        wTp(i)    =  wTp(i)/Ncount(i) 
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)

        write(3,'(16E15.7, I7)') (z_p(i)+z_p(i+1))/2.0, Ump(i), Vmp(i), Wmp(i), Tmp(i), uup(i), vvp(i), &
                                wwp(i), uvp(i), uwp(i), vwp(i), var_1(i), TTp(i),   &
                                uTp(i), vTp(i), wTp(i), Ncount(i)

      end if
    end do 
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


  if(this_proc < 2) write(*,*) 'Finished with UserCutLines_RB_conv'

  end subroutine UserCutLines_RB_conv
