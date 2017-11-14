!==============================================================================!
  subroutine Load_Ini(grid)
!------------------------------------------------------------------------------!
! This version of Load_Ini is optimised for very large meshes
! Program SUB_INI needs to be used to create files needed by this_proc
! subroutine 
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod, only: this_proc
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  real             :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer          :: j, k,  c, nearest_cell, c1, c2, s 
  integer          :: NCold  
  real,allocatable :: Xold(:),Yold(:),Zold(:)
  real,allocatable :: Uold(:),Vold(:),Wold(:),Told(:), Kold(:), Eold(:), v_2old(:), f22old(:)
  real,allocatable :: UCold(:),VCold(:),WCold(:),TCold(:), KCold(:), ECold(:), v_2Cold(:), f22Cold(:)
  real,allocatable :: UCoold(:),VCoold(:),WCoold(:),TCoold(:), KCoold(:), ECoold(:), v_2Coold(:), f22Coold(:)
  real,allocatable :: Uoold(:),Voold(:),Woold(:),Toold(:), Koold(:), Eoold(:), v_2oold(:), f22oold(:)
  real,allocatable :: UDoold(:),VDoold(:),WDoold(:),TDoold(:), KDoold(:), EDoold(:), v_2Doold(:), f22Doold(:)
  real,allocatable :: UXold(:),VXold(:),WXold(:),TXold(:), KXold(:), EXold(:), v_2Xold(:), f22Xold(:)
  real,allocatable :: UXoold(:),VXoold(:),WXoold(:),TXoold(:), KXoold(:), EXoold(:), v_2Xoold(:), f22Xoold(:)
  real             :: Us, Ws, Vs
  real             :: old_distance

  ! Variables for ReadC:
  character(len=80) :: answer
  character(len=80) :: name_in
!==============================================================================!  

  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)), '(A80)') name_in
  answer=name_in
  call To_Upper_Case(answer)
  if(answer == 'SKIP') return

  ! Save the name
  answer = name
  name = name_in

  HOTini = NO

  call Name_File(this_proc, name_in, '.ini', len_trim('.ini'))

  if(this_proc < 2) write(*,*)'now reading file:', name_in 

  open(5, file=name_in) 
  read(5,*) NCold

  allocate (Xold(NCold)); Xold = 0.0
  allocate (Yold(NCold)); Yold = 0.0
  allocate (Zold(NCold)); Zold = 0.0
  allocate (Uold(NCold)); Uold = 0.0
  allocate (Vold(NCold)); Vold = 0.0
  allocate (Wold(NCold)); Wold = 0.0
  allocate (Uoold(NCold)); Uoold = 0.0
  allocate (Voold(NCold)); Voold = 0.0
  allocate (Woold(NCold)); Woold = 0.0
  allocate (UDoold(NCold)); UDoold = 0.0
  allocate (VDoold(NCold)); VDoold = 0.0
  allocate (WDoold(NCold)); WDoold = 0.0
  allocate (UCold(NCold)); UCold = 0.0
  allocate (VCold(NCold)); VCold = 0.0
  allocate (WCold(NCold)); WCold = 0.0
  allocate (UCoold(NCold)); UCoold = 0.0
  allocate (VCoold(NCold)); VCoold = 0.0
  allocate (WCoold(NCold)); WCoold = 0.0
  allocate (UXold(NCold)); UXold = 0.0
  allocate (VXold(NCold)); VXold = 0.0
  allocate (WXold(NCold)); WXold = 0.0
  allocate (UXoold(NCold)); UXoold = 0.0
  allocate (VXoold(NCold)); VXoold = 0.0
  allocate (WXoold(NCold)); WXoold = 0.0
  if(HOTini == YES) then
    allocate (Told(NCold)); Told = 0.0
    allocate (Toold(NCold)); Toold = 0.0
    allocate (TDoold(NCold)); TDoold = 0.0
    allocate (TCold(NCold)); TCold = 0.0
    allocate (TCoold(NCold)); TCoold = 0.0
    allocate (TXold(NCold)); TXold = 0.0
    allocate (TXoold(NCold)); TXoold = 0.0
  end if
  allocate (Kold(NCold)); Kold = 0.0
  allocate (Koold(NCold)); Koold = 0.0
  allocate (KDoold(NCold)); KDoold = 0.0
  allocate (KCold(NCold)); KCold = 0.0
  allocate (KCoold(NCold)); KCoold = 0.0
  allocate (KXold(NCold)); KXold = 0.0
  allocate (KXoold(NCold)); KXoold = 0.0

  allocate (Eold(NCold)); Eold = 0.0
  allocate (Eoold(NCold)); Eoold = 0.0
  allocate (EDoold(NCold)); EDoold = 0.0
  allocate (ECold(NCold)); ECold = 0.0
  allocate (ECoold(NCold)); ECoold = 0.0
  allocate (EXold(NCold)); EXold = 0.0
  allocate (EXoold(NCold)); EXoold = 0.0

  allocate (v_2old(NCold)); v_2old = 0.0
  allocate (v_2oold(NCold)); v_2oold = 0.0
  allocate (v_2Doold(NCold)); v_2Doold = 0.0
  allocate (v_2Cold(NCold)); v_2Cold = 0.0
  allocate (v_2Coold(NCold)); v_2Coold = 0.0
  allocate (v_2Xold(NCold)); v_2Xold = 0.0
  allocate (v_2Xoold(NCold)); v_2Xoold = 0.0

  allocate (F22old(NCold)); F22old = 0.0
  allocate (F22oold(NCold)); F22oold = 0.0
  allocate (F22Doold(NCold)); F22Doold = 0.0
  allocate (F22Xold(NCold)); F22Xold = 0.0
  allocate (F22Xoold(NCold)); F22Xoold = 0.0

  j = NCold
  do k = 1, j
    if(this_proc < 2) then
      if(mod(k,20000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'  
    end if
    if(SIMULA == LES) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k)
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k)
      end if
    end if 
    if(SIMULA == ZETA.or.SIMULA == K_EPS_VV) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k)!, &
!                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
!                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
!                  v_2old(k), v_2oold(k), v_2Cold(k), v_2Coold(k), v_2Doold(k), v_2Xold(k), v_2Xoold(k), &
!                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k)
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k)!, &
!                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
!                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
!                  v_2old(k), v_2oold(k), v_2Cold(k), v_2Coold(k), v_2Doold(k), v_2Xold(k), v_2Xoold(k), &
!                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k)
      end if
    end if 
    if(SIMULA == K_EPS) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k)
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k)
      end if
    end if 
  end do
  close(5)
  if(this_proc < 2) write(*,*) 'LoaInI: finished with reading the files'

  nearest_cell = 0
  near = 0
  old_distance = HUGE
    do c = 1, grid % n_cells
      if(this_proc < 2) then
        if(mod(c,20000) == 0) write(*,*) (100.*c/(1.*grid % n_cells)), '% complete...'  
      end if
      old_distance = HUGE
      do k = 1, j
        if(Distance(Xold(k),Yold(k),Zold(k),grid % xc(c),grid % yc(c),grid % zc(c)) < old_distance) then
          old_distance = Distance(Xold(k),Yold(k),Zold(k),grid % xc(c),grid % yc(c),grid % zc(c))       
          nearest_cell =  k
        end if 
      end do  
      U % n(c)  = Uold(nearest_cell) 
      U % o(c)  = Uoold(nearest_cell) 
      U % C(c)  = UCold(nearest_cell) 
      U % Co(c) = UCoold(nearest_cell) 
      U % Do(c) = UDoold(nearest_cell) 
      U % X(c)  = UXold(nearest_cell) 
      U % Xo(c) = UXoold(nearest_cell) 
      V % n(c)  = Vold(nearest_cell) 
      V % o(c)  = Voold(nearest_cell) 
      V % C(c)  = VCold(nearest_cell) 
      V % Co(c) = VCoold(nearest_cell) 
      V % Do(c) = VDoold(nearest_cell) 
      V % X(c)  = VXold(nearest_cell) 
      V % Xo(c) = VXoold(nearest_cell) 
      W % n(c)  = Wold(nearest_cell) 
      W % o(c)  = Woold(nearest_cell) 
      W % C(c)  = WCold(nearest_cell) 
      W % Co(c) = WCoold(nearest_cell) 
      W % Do(c) = WDoold(nearest_cell) 
      W % X(c)  = WXold(nearest_cell) 
      W % Xo(c) = WXoold(nearest_cell) 
      if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==K_EPS) then
        Kin % n(c)  = Kold(near(c)) 
        Kin % o(c)  = Koold(near(c)) 
        Kin % C(c)  = KCold(near(c)) 
        Kin % Co(c) = KCoold(near(c)) 
        Kin % Do(c) = KDoold(near(c)) 
        Kin % X(c)  = KXold(near(c)) 
        Kin % Xo(c) = KXoold(near(c)) 

        Eps % n(c)  = Eold(near(c)) 
        Eps % o(c)  = Eoold(near(c)) 
        Eps % C(c)  = ECold(near(c)) 
        Eps % Co(c) = ECoold(near(c)) 
        Eps % Do(c) = EDoold(near(c)) 
        Eps % X(c)  = EXold(near(c)) 
        Eps % Xo(c) = EXoold(near(c)) 
      end if
      if(SIMULA==K_EPS_VV.or.SIMULA==ZETA) then
        v_2 % n(c)  = v_2old(near(c)) 
        v_2 % o(c)  = v_2oold(near(c)) 
        v_2 % C(c)  = v_2Cold(near(c)) 
        v_2 % Co(c) = v_2Coold(near(c)) 
        v_2 % Do(c) = v_2Doold(near(c)) 
        v_2 % X(c)  = v_2Xold(near(c)) 
        v_2 % Xo(c) = v_2Xoold(near(c)) 

        F22 % n(c)  = F22old(near(c)) 
        F22 % o(c)  = F22oold(near(c)) 
        F22 % Do(c) = F22Doold(near(c)) 
        F22 % X(c)  = F22Xold(near(c)) 
        F22 % Xo(c) = F22Xoold(near(c)) 
      end if
    end do
  do s = 1, grid % n_faces
    c1=SideC(1,s)
    c2=SideC(2,s)

    ! Interpolate density and velocity
    Us = f(s) * U % n(c1) + (1.0-f(s)) * U % n(c2)
    Vs = f(s) * V % n(c1) + (1.0-f(s)) * V % n(c2)
    Ws = f(s) * W % n(c1) + (1.0-f(s)) * W % n(c2)
    Flux(s) = (  Us * grid % sx(s)  &
               + Vs * grid % sy(s)  &
               + Ws * grid % sz(s) )
  end do 

  name = answer

  deallocate(Xold)
  deallocate(Yold)
  deallocate(Zold)
  deallocate(Uold)
  deallocate(Vold)
  deallocate(Wold)
  deallocate(Kold)
  deallocate(Eold)
  deallocate(v_2old)
  deallocate(F22old)
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(Koold)
  deallocate(Eoold)
  deallocate(v_2oold)
  deallocate(F22oold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(KDoold)
  deallocate(EDoold)
  deallocate(v_2Doold)
  deallocate(F22Doold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(KCold)
  deallocate(ECold)
  deallocate(v_2Cold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(KCoold)
  deallocate(ECoold)
  deallocate(v_2Coold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(KXold)
  deallocate(EXold)
  deallocate(v_2Xold)
  deallocate(F22Xold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
  deallocate(KXoold)
  deallocate(EXoold)
  deallocate(v_2Xoold)
  deallocate(F22Xoold)
!  deallocate(Pold)
!  deallocate(PPold)
!  deallocate(Pxold)
!  deallocate(Pyold)
!  deallocate(Pzold)
  if(HOTini==YES) then
    deallocate(Told)
    deallocate(Toold)
    deallocate(TDoold)
    deallocate(TCold)
    deallocate(TCoold)
    deallocate(TXold)
    deallocate(TXoold)
  end if

  write(*,*) 'Finished with Load_Ini  Processor: ', this_proc

  ! Restore the name
  name = answer

  end subroutine

