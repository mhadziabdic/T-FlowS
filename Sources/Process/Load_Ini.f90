!==============================================================================!
  subroutine Load_Ini(grid)
!------------------------------------------------------------------------------!
! This version of Load_Ini is optimised for very large meshes
! Program SUB_INI needs to be used to create files needed by this_proc
! subroutine 
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use allp_mod
  use Flow_Mod
  use les_mod
  use Comm_Mod, only: this_proc
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  real             :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer          :: j, k,  c, nearest_cell, c1, c2, s 
  integer          :: NCold  
  real,allocatable :: Xold(:),Yold(:),Zold(:)
  real,allocatable :: Uold(:),Vold(:),Wold(:),Told(:), Kold(:), Eold(:), v2old(:), f22old(:)
  real,allocatable :: UCold(:),VCold(:),WCold(:),TCold(:), KCold(:), ECold(:), v2Cold(:), f22Cold(:)
  real,allocatable :: UCoold(:),VCoold(:),WCoold(:),TCoold(:), KCoold(:), ECoold(:), v2Coold(:), f22Coold(:)
  real,allocatable :: Uoold(:),Voold(:),Woold(:),Toold(:), Koold(:), Eoold(:), v2oold(:), f22oold(:)
  real,allocatable :: UDoold(:),VDoold(:),WDoold(:),TDoold(:), KDoold(:), EDoold(:), v2Doold(:), f22Doold(:)
  real,allocatable :: UXold(:),VXold(:),WXold(:),TXold(:), KXold(:), EXold(:), v2Xold(:), f22Xold(:)
  real,allocatable :: UXoold(:),VXoold(:),WXoold(:),TXoold(:), KXoold(:), EXoold(:), v2Xoold(:), f22Xoold(:)
  real             :: Us, Ws, Vs
  real             :: old_distance

  ! Variables for ReadC:
  character(len=80) :: answer
  character(len=80) :: name_in
!==============================================================================!  

  call Control_Mod_Load_Initial_Solution_Name(name_in)

  answer=name_in
  call To_Upper_Case(answer)
  if(answer == 'SKIP') return

  ! Save the name
  answer = problem_name
  problem_name = name_in

  call Name_File(this_proc, name_in, '.ini')

  HOTini = NO

  if(this_proc < 2) print *,'now reading file: ', name_in 

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

  allocate (v2old(NCold)); v2old = 0.0
  allocate (v2oold(NCold)); v2oold = 0.0
  allocate (v2Doold(NCold)); v2Doold = 0.0
  allocate (v2Cold(NCold)); v2Cold = 0.0
  allocate (v2Coold(NCold)); v2Coold = 0.0
  allocate (v2Xold(NCold)); v2Xold = 0.0
  allocate (v2Xoold(NCold)); v2Xoold = 0.0

  allocate (F22old(NCold)); F22old = 0.0
  allocate (F22oold(NCold)); F22oold = 0.0
  allocate (F22Doold(NCold)); F22Doold = 0.0
  allocate (F22Xold(NCold)); F22Xold = 0.0
  allocate (F22Xoold(NCold)); F22Xoold = 0.0

  j = NCold
  do k = 1, j
    if(this_proc < 2) then
      if(mod(k,20000) == 0) print *, (100.*k/(1.*j)), '% complete...'  
    end if
    if(turbulence_model == LES) then 
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
    if(turbulence_model == K_EPS_ZETA_F .or.  &
       turbulence_model == K_EPS_V2) then 
      if(HOTini==YES) then
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k)!, &
!                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
!                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
!                  v2old(k), v2oold(k), v2Cold(k), v2Coold(k), v2Doold(k), v2Xold(k), v2Xoold(k), &
!                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k)
      else
        read(5,*) Xold(k), Yold(k), Zold(k), &
                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k)!, &
!                  Kold(k), Koold(k), KCold(k), KCoold(k), KDoold(k), KXold(k), KXoold(k), &
!                  Eold(k), Eoold(k), ECold(k), ECoold(k), EDoold(k), EXold(k), EXoold(k), &
!                  v2old(k), v2oold(k), v2Cold(k), v2Coold(k), v2Doold(k), v2Xold(k), v2Xoold(k), &
!                  F22old(k), F22oold(k), F22Doold(k), F22Xold(k), F22Xoold(k)
      end if
    end if 
    if(turbulence_model == K_EPS) then 
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
  if(this_proc < 2) print *, 'LoaInI: finished with reading the files'

  nearest_cell = 0
  old_distance = HUGE
    do c = 1, grid % n_cells
      if(this_proc < 2) then
        if(mod(c,20000) == 0) print *, (100.*c/(1.*grid % n_cells)), '% c_omplete...'  
      end if
      old_distance = HUGE
      do k = 1, j
        if(Distance(Xold(k),Yold(k),Zold(k),grid % xc(c),grid % yc(c),grid % zc(c)) < old_distance) then
          old_distance = Distance(Xold(k),Yold(k),Zold(k),grid % xc(c),grid % yc(c),grid % zc(c))       
          nearest_cell =  k
        end if 
      end do  
      U % n(c)   = Uold(nearest_cell) 
      U % o(c)   = Uoold(nearest_cell) 
      U % a(c)   = UCold(nearest_cell) 
      U % a_o(c) = UCoold(nearest_cell) 
      U % d_o(c) = UDoold(nearest_cell) 
      U % c(c)   = UXold(nearest_cell) 
      U % c_o(c) = UXoold(nearest_cell) 
      V % n(c)   = Vold(nearest_cell) 
      V % o(c)   = Voold(nearest_cell) 
      V % a(c)   = VCold(nearest_cell) 
      V % a_o(c) = VCoold(nearest_cell) 
      V % d_o(c) = VDoold(nearest_cell) 
      V % c(c)   = VXold(nearest_cell) 
      V % c_o(c) = VXoold(nearest_cell) 
      W % n(c)   = Wold(nearest_cell) 
      W % o(c)   = Woold(nearest_cell) 
      W % a(c)   = WCold(nearest_cell) 
      W % a_o(c) = WCoold(nearest_cell) 
      W % d_o(c) = WDoold(nearest_cell) 
      W % c(c)   = WXold(nearest_cell) 
      W % c_o(c) = WXoold(nearest_cell) 
      if(turbulence_model == K_EPS_V2 .or.  &
         turbulence_model == K_EPS_ZETA_F     .or.  &
         turbulence_model == K_EPS) then
        Kin % n(c)   = Kold(nearest_cell) 
        Kin % o(c)   = Koold(nearest_cell) 
        Kin % a(c)   = KCold(nearest_cell) 
        Kin % a_o(c) = KCoold(nearest_cell) 
        Kin % d_o(c) = KDoold(nearest_cell) 
        Kin % c(c)   = KXold(nearest_cell) 
        Kin % c_o(c) = KXoold(nearest_cell) 

        Eps % n(c)   = Eold(nearest_cell) 
        Eps % o(c)   = Eoold(nearest_cell) 
        Eps % a(c)   = ECold(nearest_cell) 
        Eps % a_o(c) = ECoold(nearest_cell) 
        Eps % d_o(c) = EDoold(nearest_cell) 
        Eps % c(c)   = EXold(nearest_cell) 
        Eps % c_o(c) = EXoold(nearest_cell) 
      end if
      if(turbulence_model == K_EPS_V2 .or.  &
         turbulence_model == K_EPS_ZETA_F) then
        v2 % n(c)   = v2old(nearest_cell) 
        v2 % o(c)   = v2oold(nearest_cell) 
        v2 % a(c)   = v2Cold(nearest_cell) 
        v2 % a_o(c) = v2Coold(nearest_cell) 
        v2 % d_o(c) = v2Doold(nearest_cell) 
        v2 % c(c)   = v2Xold(nearest_cell) 
        v2 % c_o(c) = v2Xoold(nearest_cell) 

        F22 % n(c)   = F22old(nearest_cell) 
        F22 % o(c)   = F22oold(nearest_cell) 
        F22 % d_o(c) = F22Doold(nearest_cell) 
        F22 % c(c)   = F22Xold(nearest_cell) 
        F22 % c_o(c) = F22Xoold(nearest_cell) 
      end if
    end do
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Interpolate density and velocity
    Us = grid % f(s) * U % n(c1) + (1.0 - grid % f(s)) * U % n(c2)
    Vs = grid % f(s) * V % n(c1) + (1.0 - grid % f(s)) * V % n(c2)
    Ws = grid % f(s) * W % n(c1) + (1.0 - grid % f(s)) * W % n(c2)
    flux(s) = (  Us * grid % sx(s)  &
               + Vs * grid % sy(s)  &
               + Ws * grid % sz(s) )
  end do 

  problem_name = answer

  deallocate(Xold)
  deallocate(Yold)
  deallocate(Zold)
  deallocate(Uold)
  deallocate(Vold)
  deallocate(Wold)
  deallocate(Kold)
  deallocate(Eold)
  deallocate(v2old)
  deallocate(F22old)
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(Koold)
  deallocate(Eoold)
  deallocate(v2oold)
  deallocate(F22oold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(KDoold)
  deallocate(EDoold)
  deallocate(v2Doold)
  deallocate(F22Doold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(KCold)
  deallocate(ECold)
  deallocate(v2Cold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(KCoold)
  deallocate(ECoold)
  deallocate(v2Coold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(KXold)
  deallocate(EXold)
  deallocate(v2Xold)
  deallocate(F22Xold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
  deallocate(KXoold)
  deallocate(EXoold)
  deallocate(v2Xoold)
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

  print *, 'Finished with Load_Ini  Processor: ', this_proc

  ! Restore the name
  problem_name = answer

  end subroutine

