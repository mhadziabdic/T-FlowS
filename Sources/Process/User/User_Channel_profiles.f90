!======================================================================!
  subroutine User_Channel_profiles(grid, y) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous directions. 
! The results are writen in files name_res.dat and name_res_plus.dat
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
!  real :: y(-NbC:NC)
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
    integer             :: Nprob, pl, c, i, count, s, c1, c2
    character           :: namCoo*80, namPro*80, namRes*80
    character           :: namRes_plus*80
    real,allocatable    :: z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:), ind(:),       & 
                                 var_1(:), var_2(:), var_3(:), &
                                 Wall_p(:), Ufric_p(:)
    integer,allocatable :: Np(:), Ncount(:)
    real                :: Lscale, Twall, Tfric, Dwall
    logical             :: there
!======================================================================!

    namPro = name

    namRes = name
    namRes(len_trim(name)+1:len_trim(name)+8) = '_res.dat'
    namRes_plus = name
    namRes_plus(len_trim(name)+1:len_trim(name)+13) = '_res_plus.dat'
!>>>>>>>>>>>>>>>>>>>>>>!
!     read 1D file     !
!>>>>>>>>>>>>>>>>>>>>>>!
    namCoo = name
    namCoo(len_trim(name)+1:len_trim(name)+3) = '.1D'

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

    if(this_proc < 2)  write(6, *) '# Now reading the file:', namCoo
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
    allocate(var_1(Nprob)); var_1=0.0
    allocate(var_2(Nprob)); var_2=0.0
    allocate(var_3(Nprob)); var_3=0.0
    allocate(Ufric_p(Nprob)); Ufric_p = 0.0

    allocate(Ncount(Nprob)); Ncount=0
    count = 0
    if(HOT==YES) then
      allocate(Tmp(Nprob));   Tmp=0.0
      allocate(TTp(Nprob));   TTp=0.0
      allocate(uTp(Nprob));   uTp=0.0
      allocate(vTp(Nprob));   vTp=0.0
      allocate(wTp(Nprob));   wTp=0.0
    end if  

    Lscale = 0.0
!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
    do i = 1, Nprob-1
      do c=1, grid % n_cells
        Lscale = max(Lscale,WallDs(c))
        if(y(c) > (z_p(i+1)) .and. y(c) < (z_p(i))) then
          if(ROT == YES) then
            Wall_p(i) = Wall_p(i) + y(c) 
          else
            Wall_p(i) = Wall_p(i) + WallDs(c) 
          end if 
          if(SIMULA==LES.or.SIMULA==DNS.or.SIMULA==HYB_ZETA.or.SIMULA==DES_SPA) then       
            Ump(i)   = Ump(i) + U % mean(c)
            Vmp(i)   = Vmp(i) + V % mean(c)
            Wmp(i)   = Wmp(i) + W % mean(c)
          else
            Ump(i)   = Ump(i) + U % n(c)
            Vmp(i)   = Vmp(i) + V % n(c)
            Wmp(i)   = Wmp(i) + W % n(c)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==K_EPS.or.SIMULA==HYB_ZETA) then
            if(IsNearWall(c)) then
              if(Ynd(c) > 5.0) then        
                Ufric_p(i) = Ufric_p(i) + &
                sqrt(max(abs(bulk(1) % p_drop_x), abs(bulk(1) % p_drop_y),&
                abs(bulk(1) % p_drop_z))/DENc(material(c))) 
              else  
                Ufric_p(i) = Ufric_p(i) + sqrt((VISc*(U % n(c)**2 + V % n(c)**2 + &
                           W% n(c)**2)**0.5/WallDs(c))/DENc(material(c)))
              end if 
            end if  
            uup(i)   = uup(i) + Kin%n(c)
            vvp(i)   = vvp(i) + Eps%n(c)
            wwp(i)   = wwp(i) + VISt(c)*(U % z(c) + W % x(c)) 
            uvp(i)   = uvp(i) + VISt(c)/VISc
          end if
          if(SIMULA==ZETA.or.SIMULA==K_EPS_VV) then
            uwp(i)   = uwp(i) + f22%n(c)
            vwp(i)   = vwp(i) + v_2%n(c) 
          end if
          if(SIMULA==EBM.or.SIMULA==HJ) then
            uup(i)   = uup(i) + uu%n(c) 
            vvp(i)   = vvp(i) + vv%n(c) 
            wwp(i)   = wwp(i) + ww%n(c)
            uvp(i)   = uvp(i) + uv%n(c)
            uwp(i)   = uwp(i) + uw%n(c)
            vwp(i)   = vwp(i) + vw%n(c)
            var_1(i) = var_1(i) + Kin%n(c)
            var_2(i) = var_2(i) + Eps%n(c)
            if(SIMULA==EBM) var_3(i) = var_3(i) + f22%n(c) 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + sqrt((VISc*sqrt(U % n(c)**2 +&
              V % n(c)**2 + W % n(c)**2)/WallDs(c))/DENc(material(c)))
            end if
          end if 
          if(SIMULA==SPA_ALL) then
            uup(i)   = uup(i) + 2.0*VISt(c)*U % x(c) 
            vvp(i)   = vvp(i) + 2.0*VISt(c)*V % y(c)
            wwp(i)   = wwp(i) + 2.0*VISt(c)*W % z(c)
            uvp(i)   = uvp(i) + VISt(c)*(U % y(c) + V % x(c))
            uwp(i)   = uwp(i) + VISt(c)*(U % z(c) + W % x(c))
            vwp(i)   = vwp(i) + VISt(c)*(V % z(c) + W % y(c))
            var_1(i) = var_1(i) + VISt(c)/VISc
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + sqrt((VISc*sqrt(U % n(c)**2 +&
              V % n(c)**2 + W % n(c)**2)/WallDs(c))/DENc(material(c)))
            end if 
          end if
          if(SIMULA==LES.or.SIMULA==DES_SPA.or.SIMULA == DNS.or.SIMULA==HYB_ZETA) then
            uup(i)   = uup(i) + (uu % mean(c)- U % mean(c) * U % mean(c))
            vvp(i)   = vvp(i) + (vv % mean(c)- V % mean(c) * V % mean(c))
            wwp(i)   = wwp(i) + (ww % mean(c)- W % mean(c) * W % mean(c))
            uvp(i)   = uvp(i) + (uv % mean(c)- U % mean(c) * V % mean(c))
            uwp(i)   = uwp(i) + (uw % mean(c)- U % mean(c) * W % mean(c))
            vwp(i)   = vwp(i) + (vw % mean(c)- V % mean(c) * W % mean(c))
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + sqrt((VISc*sqrt(U % mean(c)**2 +&
              V % mean(c)**2 + W % mean(c)**2)/WallDs(c))/DENc(material(c)))
            end if 
          end if
          if(SIMULA==HYB_ZETA) then
            var_1(i) = var_1(i) + Kin % n(c) 
            var_2(i) = var_2(i) + VISt(c)/VISc 
            var_3(i) = var_3(i) + VISt_sgs(c)/VISc 
            if(IsNearWall(c)) then
              Ufric_p(i) = Ufric_p(i) + sqrt((VISc*sqrt(U % mean(c)**2 +&
              V % mean(c)**2 + W % mean(c)**2)/WallDs(c))/DENc(material(c)))
            end if 
          end if 
 
          if(HOT==YES) then
            if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS &
               .or.SIMULA==HYB_PITM.or.SIMULA==HYB_ZETA) then
              Tmp(i)   = Tmp(i) + T % mean(c)
              TTp(i)   = TTp(i) + (TT % mean(c) - T % mean(c) * T % mean(c))
              uTp(i)   = uTp(i) + (uT % mean(c) - u % mean(c) * T % mean(c))
              vTp(i)   = vTp(i) + (vT % mean(c) - v % mean(c) * T % mean(c))
              wTp(i)   = wTp(i) + (wT % mean(c) - w % mean(c) * T % mean(c))
            else
              Tmp(i)   = Tmp(i) + T % n(c)
              uTp(i)   = uTp(i) + uT % n(c) 
              vTp(i)   = vTp(i) + vT % n(c)
              wTp(i)   = wTp(i) + wT % n(c)
            end if
          end if
          Ncount(i) = Ncount(i) + 1
        end if
      end do 
    end do 

   if(HOT==YES) then
     Dwall = 0.0
     do c = 1, grid % n_cells
       if(WallDs(c) > Dwall) then
         Dwall = WallDs(c)
         Tinf  = T % n(c)
       end if
     end do

     call wait

     if(Qflux> 0.0) then
       call GloMin(Tinf)
     else
       call GloMax(Tinf)
     end if
    
     do s = 1, grid % n_faces
       c1=grid % faces_c(1,s)
       c2=grid % faces_c(2,s)
       if(c2  < 0) then
         if( TypeBC(c2) ==  WALL.or. TypeBC(c2) ==  WALLFL) then
           Twall = T % n(c2) 
           Numax = T % q(c2)/(CONc(material(c2))*(Twall-Tinf))
         end if
       end if
     end do
   end if

!---- average over all processors
    do pl=1, Nprob-1
      call IGlSum(Ncount(pl))

      call GloSum(Wall_p(pl))

      call GloSum(Ump(pl))
      call GloSum(Vmp(pl))
      call GloSum(Wmp(pl))

      call GloSum(var_1(pl))
      call GloSum(var_2(pl))
      call GloSum(var_3(pl))

      call GloSum(uup(pl))
      call GloSum(vvp(pl))
      call GloSum(wwp(pl))

      call GloSum(uvp(pl))
      call GloSum(uwp(pl))
      call GloSum(vwp(pl))
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
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        Ufric_p(i)=  Ufric_p(i)/Ncount(i)
        if(HOT == YES) then
          Tmp(i)    =  Tmp(i)/Ncount(i)
          TTp(i)    =  TTp(i)/Ncount(i)
          uTp(i)    =  uTp(i)/Ncount(i)
          vTp(i)    =  vTp(i)/Ncount(i)
          wTp(i)    =  wTp(i)/Ncount(i)
        end if
      end if
    end do 


    if(SIMULA == LES.or.SIMULA==DES_SPA.or.SIMULA==DNS) then
      do i = 1, (Nprob-1)/2
        Ump(i)   =  (Ump(i)+Ump(Nprob-i))/2.0
        Vmp(i)   =  (Vmp(i)+Vmp(Nprob-i))/2.0 
        Wmp(i)   =  (Wmp(i)+Wmp(Nprob-i))/2.0
        uup(i)   =  (uup(i)+uup(Nprob-i))/2.0
        vvp(i)   =  (vvp(i)+vvp(Nprob-i))/2.0
        wwp(i)   =  (wwp(i)+wwp(Nprob-i))/2.0
        uvp(i)   =  (uvp(i)+abs(uvp(Nprob-i)))/2.0
        uwp(i)   =  (uwp(i)+abs(uwp(Nprob-i)))/2.0
        vwp(i)   =  (vwp(i)+vwp(Nprob-i))/2.0
        Ufric_p(i) = (Ufric_p(i) + Ufric_p(Nprob-i))/2
      end do
    end if

    Ufric = 0.0
    do i = 1, Nprob
      Ufric = max(Ufric_p(i),Ufric)
    end do
    Tfric = Qflux/(DENc(1)*CAPc(1)*Ufric)


    if(Ufric == 0.0) then
      if(this_proc < 2) write(*,*) 'Friction velocity is zero in UserCutLines_channel.f90 !'
      return
    end if

    open(3,file=namRes)
    open(4,file=namRes_plus)

    do i = 3, 4
    write(i,'(A1,2(A8,F12.5, 3X))') '#', 'Utau = ', &
    Ufric, 'Re_tau = ', Ufric/VISc
    if(HOT==YES) write(i,'(A1,3(A10, F10.5, 2X))') '#', 'Qflux = ', abs(Qflux), &
    'K+ = q/Uf =', abs(Qflux)/Ufric*(VISc/CONc(material(1)))**ONE_THIRD,'Nu max = ', &
    Numax
    if(SIMULA == DNS.or.SIMULA==LES.or.SIMULA==HYB_ZETA.or.&
       SIMULA==DES_SPA) then
      if(HOT == YES) then 
        write(i,'(A1,2X,A100)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, &
                                 6:uu, 7:vv, 8:ww, '//  &
                  '9:uv, 10:uw, 11:vw, 12: tt, 13: uT, &
                   14: vT, 15: wT, 16:Kin, 17: VISt/VISc' 
      else
        write(i,'(A1,2X,A64)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:uu, 6:vv, &
        7:ww, 8:uv, 9:uw, 10:vw, 11:Kin 12: VISt/VISc' 
      end if
    else if(SIMULA==K_EPS) then
      if(HOT == YES) then 
        write(i,'(A1,2X,A100)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:Kin 7:Eps, 8:uw, '//  &
                  '9:VISt/VISc, 10: uT, 11: vT, 12: wT, 13: VISt/VISc' 
      else
        write(i,'(A1,2X,A70)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, 7:uw, 8:VISt/VISc' 
      end if
    else if(SIMULA==ZETA.or.SIMULA==K_EPS_VV) then
      if(HOT == YES) then 
        write(i,'(A1,2X,A100)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:Kin, 7:Eps, 8:uw, '//  &
                  '9:VISt/VISc, 10:f22, 11:v2, 12:uT, 13:vT, 14:wT ' 
      else
        write(i,'(A1,2X,A64)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:Kin, 6:Eps, &
        7:uw, 8:VISt/VISc, 9:f22, 10:v2' 
      end if
    else if(SIMULA==HJ.or.SIMULA==EBM) then
      if(HOT == YES) then 
        write(i,'(A1,2X,A100)') '#', '1:Xrad, 2:U, 3:V, 4:W, 5:T, 6:uu, 7:vv, 8:ww, '//  &
                  '9:uv, 10:uw, 11:vw, 12: uT, 13: vT, 45: wT, 15:Kin, 16: VISt/VISc' 
      else
        write(i,'(A1,2X,A64)') '#', '1:Xrad, 2:U, 3:V, 4:W, &
        5:uu, 6:vv, 7:ww, 8:uv, 9:uw, 10:vw, 11:Kin' 
      end if
    end if
    end do ! end i

    if(HOT==YES) then
      if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS &
      .or.SIMULA==HYB_PITM.or.SIMULA==HYB_ZETA) then
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(18E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), TTp(i), uTp(i), vTp(i), wTp(i), &
            var_1(i), var_2(i), var_3(i)
          end if
        end do 
      else
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(3,'(17E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), uTp(i), vTp(i), wTp(i), &
            var_1(i), var_2(i), var_3(i)
          end if
        end do 
      end if
    else 
      do i = 1, Nprob
        if(Ncount(i) /= 0) then
          write(3,'(13E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), var_1(i), var_2(i), var_3(i)
        end if
      end do 
    end if
    close(3)

    do i = 1, Nprob-1
      Wall_p(i)= Wall_p(i)*Ufric/VISc 
      Ump(i)   = Ump(i)/Ufric 
      Vmp(i)   = Vmp(i)/Ufric 
      Wmp(i)   = Wmp(i)/Ufric 
      if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==K_EPS.or.SIMULA==HYB_ZETA) then
        uup(i)   = uup(i)/Ufric**2          ! Kin%n(c)
        vvp(i)   = vvp(i)*VISc/Ufric**4.0   ! Eps%n(c)
        wwp(i)   = wwp(i)/Ufric**2          ! VISt(c)*(U % z(c) + W % x(c)) 
        uvp(i)   = uvp(i)                   ! VISt(c)/VISc
      else if(SIMULA==ZETA.or.SIMULA==K_EPS_VV) then
        uwp(i)   = uwp(i)*VISc/Ufric**2.0   ! f22%n(c)
        vwp(i)   = vwp(i)                   ! v_2%n(c) 
      else if(SIMULA==EBM.or.SIMULA==HJ) then
        uup(i)   = uup(i)/Ufric**2          ! uu%n(c) 
        vvp(i)   = vvp(i)/Ufric**2          ! vv%n(c) 
        wwp(i)   = wwp(i)/Ufric**2          ! ww%n(c)
        uvp(i)   = uvp(i)/Ufric**2          ! uv%n(c)
        uwp(i)   = uwp(i)/Ufric**2          ! uw%n(c)
        vwp(i)   = vwp(i)/Ufric**2          ! vw%n(c)
        var_1(i) = var_1(i)/Ufric**2        ! Kin%n(c)
        var_2(i) = var_2(i)*VISc/Ufric**4.0 ! Eps%n(c)
      else if(SIMULA==SPA_ALL) then
        uup(i)   = uup(i)/Ufric**2   ! 2.0*VISt(c)*U % x(c) 
        vvp(i)   = vvp(i)/Ufric**2   ! 2.0*VISt(c)*V % y(c)
        wwp(i)   = wwp(i)/Ufric**2   ! 2.0*VISt(c)*W % z(c)
        uvp(i)   = uvp(i)/Ufric**2   ! VISt(c)*(U % y(c) + V % x(c))
        uwp(i)   = uwp(i)/Ufric**2   ! VISt(c)*(U % z(c) + W % x(c))
        vwp(i)   = vwp(i)/Ufric**2   ! VISt(c)*(V % z(c) + W % y(c))
        var_1(i) = var_1(i)          ! VISt(c)/VISc
      else if(SIMULA==LES.or.SIMULA==DES_SPA.or.&
              SIMULA == DNS.or.SIMULA==HYB_ZETA) then
        uup(i)   = uup(i)/Ufric**2 ! (uu%mean(c)-U%mean(c)*U%mean(c))
        vvp(i)   = vvp(i)/Ufric**2 ! (vv%mean(c)-V%mean(c)*V%mean(c))
        wwp(i)   = wwp(i)/Ufric**2 ! (ww%mean(c)-W%mean(c)*W%mean(c))
        uvp(i)   = uvp(i)/Ufric**2 ! (uv%mean(c)-U%mean(c)*V%mean(c))
        uwp(i)   = uwp(i)/Ufric**2 ! (uw%mean(c)-U%mean(c)*W%mean(c))
        vwp(i)   = vwp(i)/Ufric**2 ! (vw%mean(c)-V%mean(c)*W%mean(c))
      else if(SIMULA==HYB_ZETA) then
        var_1(i) = var_1(i)/Ufric**2        ! Kin % n(c) 
        var_2(i) = var_2(i)                 ! VISt(c)/VISc 
        var_3(i) = var_3(i)                 ! VISt_sgs(c)/VISc 
      end if 
 
      if(HOT==YES) then
        if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS &
          .or.SIMULA==HYB_PITM.or.SIMULA==HYB_ZETA) then
          Tmp(i)   = (Twall - Tmp(i))/Tfric    ! T % mean(c)
          TTp(i)   = TTp(i)/Tfric**2     !(TT%mean(c)-T%mean(c)*T%mean(c))
          uTp(i)   = uTp(i)/(Ufric*Tfric)!(uT%mean(c)-u%mean(c)*T%mean(c))
          vTp(i)   = vTp(i)/(Ufric*Tfric)!(vT%mean(c)-v%mean(c)*T%mean(c))
          wTp(i)   = wTp(i)/(Ufric*Tfric)!(wT%mean(c)-w%mean(c)*T%mean(c))
        else
          Tmp(i)   = (Twall - Tmp(i))/Tfric ! T % n(c)
          uTp(i)   = uTp(i)/(Ufric*Tfric)   ! uT % n(c) 
          vTp(i)   = vTp(i)/(Ufric*Tfric)   ! vT % n(c)
          wTp(i)   = wTp(i)/(Ufric*Tfric)   ! wT % n(c)
        end if
      end if
    end do 


    if(HOT==YES) then
      if(SIMULA == LES.or.SIMULA == DES_SPA.or.SIMULA == DNS &
      .or.SIMULA==HYB_PITM.or.SIMULA==HYB_ZETA) then
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(4,'(18E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), TTp(i), uTp(i), vTp(i), wTp(i), &
            var_1(i), var_2(i), var_3(i)
          end if
        end do 
      else
        do i = 1, Nprob
          if(Ncount(i) /= 0) then
            write(4,'(17E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i), Tmp(i), &
            uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
            vwp(i), uTp(i), vTp(i), wTp(i), &
            var_1(i), var_2(i), var_3(i)
          end if
        end do 
      end if
    else 
      do i = 1, Nprob
        if(Ncount(i) /= 0) then
          write(4,'(13E15.7)') Wall_p(i), Ump(i), Vmp(i), Wmp(i),   &
          uup(i), vvp(i), wwp(i), uvp(i), uwp(i),      &
          vwp(i), var_1(i), var_2(i), var_3(i)
        end if
      end do 
    end if

    close(4)

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
    if(HOT==YES) then
      deallocate(Tmp)
      deallocate(TTp)
      deallocate(uTp)
      deallocate(vTp)
      deallocate(wTp)
    end if
  end subroutine User_Channel_profiles
