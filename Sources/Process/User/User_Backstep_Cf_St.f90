!======================================================================!
  subroutine User_Backstep_Cf_St(grid) 
!----------------------------------------------------------------------!
! Subroutine extracts skin friction coefficient and Stanton number for !
! backstep case.                                                       !
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
  use Grid_Mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
!  real :: y(-NbC:NC)
  type(Grid_Type) :: grid
  real :: Rad_2, Ufric 
!------------------------------[Calling]-------------------------------!
  integer             :: Nprob, pl, c, dummy, i, count, k, c1, c2, s, N_hor
  character           :: namCoo*80, namPro*80, answer*80, JetIn*17, namOut*80
  real,allocatable    :: r1_p(:), r2_p(:), z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 uum(:), vvm(:), wwm(:), &
                                 uvm(:), uwm(:), vwm(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 Var_1(:), Var_2(:), Var_3(:), Rad_mp(:), &
                                 Var_4(:), Var_5(:)  
  integer,allocatable :: Np(:), Ncount(:)
  real                :: R, Urad_n, Utan_n, R1, R2, Urad, Utan 
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    if(this_proc < 2) write(6, *) '# Now reading the file: X_coordinate.dat ' 
    open(9, file='x_coordinate.dat')

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
    allocate(Var_1(Nprob));  Var_1=0.0
    allocate(Var_2(Nprob));  Var_2=0.0
    allocate(Var_3(Nprob));  Var_3=0.0
    allocate(Var_4(Nprob));  Var_4=0.0
    allocate(Var_5(Nprob));  Var_5=0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0

    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do i = 1, Nprob-1
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2).eq.WALLFL.and.T % q(c2) > 1.0e-8) then
          if(grid % xc(c1) > z_p(i) .and. grid % xc(c1) < z_p(i+1)) then
            Ump(i)   = Ump(i) + U % n(c1)
            Vmp(i)   = Vmp(i) + V % n(c1)
            Wmp(i)   = Wmp(i) + W % n(c1)
            Var_1(i) = Var_1(i) + 2.0*(VISc*U%n(c1)/WallDs(c1))/11.3**2 
            Var_2(i) = Var_2(i) + sqrt(abs(U%n(c1))*WallDs(c1)/VISc) 
            Var_3(i) = Var_3(i) + 0.1/((T%n(c2)-20)*11.3) 
            Var_5(i) = Var_5(i) + T % n(c2) 
            Rad_mp(i) = Rad_mp(i) + grid % zc(c1) 
            Ncount(i) = Ncount(i) + 1
          end if
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
    call GloSum(Var_1(pl))
    call GloSum(Var_2(pl))
    call GloSum(Var_3(pl))
    call GloSum(Var_4(pl))
    call GloSum(Var_5(pl))

    count =  count + Ncount(pl) 

    if(HOT==YES) then
      call GloSum(Tmp(pl))
      call GloSum(TTp(pl))
      call GloSum(uTp(pl))
      call GloSum(vTp(pl))
      call GloSum(wTp(pl))
    end if
  end do

  JetIn  = 'Result_Cf_St'
  write(JetIn(13:17),'(A4)') '.dat'

  open(3,file=JetIn)
  write(3,*) '# x, Cf, St, U, T, yPlus'
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
      Var_1(i)  =  Var_1(i)/Ncount(i)
      Var_2(i)  =  Var_2(i)/Ncount(i)
      Var_3(i)  =  Var_3(i)/Ncount(i)
      Var_4(i)  =  Var_4(i)/Ncount(i)
      Var_5(i)  =  Var_5(i)/Ncount(i)
      Rad_mp(i) =  Rad_mp(i)/Ncount(i)

      write(3,'(6E15.7)') (z_p(i)+z_p(i+1))/(2.*0.038), &
      Var_1(i), Var_3(i), Ump(i), Var_5(i), var_2(i) 

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
  deallocate(uum)
  deallocate(vvm)
  deallocate(wwm)
  deallocate(uvm)
  deallocate(uwm)
  deallocate(vwm)
  deallocate(Var_1)
  deallocate(Var_2)
  deallocate(Var_3)
  deallocate(Var_4)
  deallocate(Var_5)
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


  if(this_proc < 2) write(*,*) 'Finished with User_Backstep_Cf_St'

  end subroutine
