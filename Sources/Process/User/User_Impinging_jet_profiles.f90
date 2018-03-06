!======================================================================!
  SUBROUTINE User_Impinging_jet_profiles(grid, y) 
!----------------------------------------------------------------------!
! Subroutine reads ".1D" file created by the "Generator" or Neu2TFlowS !
! and extracts profiles on several locations that corresponds with the !
! experimental measurements.                                           !
!----------------------------------------------------------------------!
    use all_mod
    use allp_mod
    use les_mod
    use pro_mod
    use par_mod
    use rans_mod
    use Grid_Mod
    use Parameters_Mod
!----------------------------------------------------------------------!
    implicit none
!-----------------------------[Arguments]------------------------------!
    type(Grid_Type) :: grid
!-----------------------------[Arguments]------------------------------!
    real :: y(-grid % n_bnd_cells:grid % n_cells)
    real :: Ufric
!------------------------------[Calling]-------------------------------!
    interface
      logical function Approx(A,B,tol)
        real           :: A,B
        real, optional :: tol
      end function Approx
    end interface
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: Nprob, pl, c, dummy, i, count, k, c1, c2, s
  CHARACTER           :: namCoo*80, namPro*80, answer*80, JetIn*16, namOut*80
  REAL,ALLOCATABLE    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_mp(:), &
                                 var_4(:), var_5(:)  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, Urad_mean, Utan_mean, R1, R2, Urad, Utan, lnum
  REAL                :: Rad_2, Uaver
  logical             :: there
!--------------------------------[CVS]---------------------------------!
!  $Id: UserProbe1D.f90,v 1.16 2002/11/25 10:33:17 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserProbe1D.f90,v $  
!======================================================================!
!    pi = 3.141592
  Uaver = 1.14

  namCoo = name
  namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'
  if(this_proc < 2) write(6, *) '# Now reading the file:', namCoo

  INQUIRE( file=namCoo, EXIST=THERE ) 
  if(.NOT.THERE) then
    if(this_proc < 2) write(*,*) "==================================================================="
    if(this_proc < 2) write(*,*) "In order to extract profiles and write them in ascii files"
    if(this_proc < 2) write(*,*) "the code has to read cell-faces coordinates "
    if(this_proc < 2) write(*,*) "in wall-normal direction in the ascii file 'case_name'.1D."
    if(this_proc < 2) write(*,*) "The file format should be as follows:"
    if(this_proc < 2) write(*,*) "10  ! number of cells + 1"
    if(this_proc < 2) write(*,*) "1 0.0"
    if(this_proc < 2) write(*,*) "2 0.1"
    if(this_proc < 2) write(*,*) "3 0.2"
    if(this_proc < 2) write(*,*) "... "
    if(this_proc < 2) write(*,*) "==================================================================="
    return
  end if

  open(9, FILE=namCoo)

!---- write the number of searching intervals
  read(9,*) Nprob
  allocate(z_p(Nprob*2))
  allocate(ind(Nprob*2))
!---- read the intervals positions
  do pl=1,Nprob
    read(9,*) ind(pl), z_p(pl)
  end do
  close(9)
  call SSORT (z_p, ind, Nprob,0)

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
  allocate(Rad_mp(Nprob));  Rad_mp=0.0

  allocate(Ncount(Nprob)); Ncount=0
  count = 0

  if(HOT==YES) then
    allocate(Tmp(Nprob));   Tmp=0.0
  end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!

  do k = 0, 6
    if(k == 0) then
      R1 = 0.04   
      R2 = 0.0
      lnum = 0.0
    else if(k == 1) then
      R1 = 1.0    
      R2 = 0.992  
      lnum = 0.5
    else if(k == 2) then
      R1 = 2.1500 
      R2 = 2.0 
      lnum = 1.0
    else if(k == 3) then
      R1 = 3.0684
      R2 = 2.9744
      lnum = 1.5
    else if(k == 4) then
      R1 = 4.1433 
      R2 = 3.9098
      lnum = 2.0
    else if(k == 5) then
      R1 = 0.5347000E+01 
      R2 = 0.4803200E+01 
      lnum = 2.5
    else if(k == 6) then
      R1 = 0.6000000E+01
      R2 = 0.5876600E+01
      lnum = 3.0
    end if  

    do i = 1, Nprob-1
      do c = 1, grid % n_cells
        Rad_2  = (grid % xc(c)**2 + grid % yc(c)**2)**0.5 + tiny
        if(Rad_2 < R1 .and. Rad_2 > R2) then
          if(grid % zc(c) < z_p(i).and.grid % zc(c) > z_p(i+1)) then
            R   = (grid % xc(c)**2 + grid % yc(c)**2)**0.5 + tiny
            Urad   = (U % n(c)*grid % xc(c)/R + &
                      V % n(c)*grid % yc(c)/R)
            Utan   = (-U % n(c)*grid % yc(c)/R  + &
                       V % n(c)*grid % xc(c)/R) 
            Ump(i)   = Vmp(i) + (U % n(c)**2.0 + V % n(c)** 2.0 &
                     + W % n(c)**2.0)**0.5 
            Vmp(i)   = Ump(i) + Urad
            Wmp(i)   = Wmp(i) + W % n(c)

            if(SIMULA == K_EPS.or.SIMULA==ZETA) then
              uuP(i)   = uup(i) + Kin % n(c)  
              vvp(i)   = vvp(i) + Eps % n(c)
              vwp(i)   = vwp(i) + VISt(c)/VISc
            end if
            if(SIMULA == ZETA) then  
              Wmp(i)   = Wmp(i) + W % n(c)
              wwp(i)   = wwp(i) + v_2 % n(c)
              uvp(i)   = uvp(i) + f22 % n(c)
            end if
            if(HOT==YES) then
              Tmp(i)   = Tmp(i) + T % n(c)
            end if
     
            Rad_mp(i) = Rad_mp(i) + R
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

      count =  count + Ncount(pl) 

      if(HOT==YES) then
        call GloSum(Tmp(pl))
      end if
    end do

    JetIn  = 'Jet_res_x'
    write(JetIn(9:16),'(F3.1,A4)') lnum, '.dat'

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Wmp(i)    =  Wmp(i)/Ncount(i)
        Ump(i)    =  Ump(i)/Ncount(i)
        Vmp(i)    =  Vmp(i)/Ncount(i)
        uup(i)    =  uup(i)/Ncount(i)
        vvp(i)    =  vvp(i)/Ncount(i)
        wwp(i)    =  wwp(i)/Ncount(i)
        uvp(i)    =  uvp(i)/Ncount(i)
        Tmp(i)    =  Tmp(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
      end if
    end do

    open(3,FILE=JetIn)
    write(3,'(A1,2X,A100)') '#', '1:Xrad, 2:Umag, 3:Urad, 4:Uaxi, 5:Kin, &
                                  6:Eps, 7:Temp, 8:VISt/VISc, '//  &
                                 '9:v_2, 10:f22'  

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        write(3,'(9E15.7)') (z_p(i)+z_p(i+1))/4.0, Ump(i)/Uaver, &
        Vmp(i)/Uaver, Wmp(i)/Uaver, uup(i)/Uaver**2, vvp(i), Tmp(i), vwp(i), &
        wwp(i), uvp(i) 
      end if
    end do 
    close(3)

    do i = 1, Nprob
      Ncount(i) = 0
      Wmp(i)    = 0.0 
      Ump(i)    = 0.0 
      Vmp(i)    = 0.0 
      uup(i)    = 0.0 
      vvp(i)    = 0.0 
      wwp(i)    = 0.0 
      uvp(i)    = 0.0 
      Tmp(i)    = 0.0 
      Rad_mp(i) =  0.0
    end do
    if(this_proc < 2) write(*,*) 'Finished with profile r/D =  ', lnum
  end do   !end number of radius

  deallocate(Np)
  deallocate(z_p)
  deallocate(Ump)
  deallocate(Vmp)
  deallocate(Wmp)
  deallocate(uup)
  deallocate(vvp)
  deallocate(wwp)
  deallocate(uvp)
  deallocate(Rad_mp)
  deallocate(Ncount)
  if(HOT==YES) then
    deallocate(Tmp)
  end if

  if(this_proc < 2) write(*,*) 'Finished with User_Impinging_jet_profiles'

  END SUBROUTINE User_Impinging_jet_profiles
