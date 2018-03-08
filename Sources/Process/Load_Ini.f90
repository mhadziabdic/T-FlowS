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
  use Tokenizer_Mod
  use Grid_Mod
  use Parameters_Mod
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

  call Tokenizer_Mod_Read_Line(CMN_FILE)
  read(line % tokens(1), '(A80)') name_in
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

  j = NCold
  do k = 1, j
    if(this_proc < 2) then
      if(mod(k,20000) == 0) write(*,*) (100.*k/(1.*j)), '% c_omplete...'  
    end if
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
  end do
  close(5)
  if(this_proc < 2) write(*,*) 'LoaInI: finished with reading the files'

  nearest_cell = 0
  near = 0
  old_distance = HUGE
  do c = 1, grid % n_cells
    if(this_proc < 2) then
      if(mod(c,20000) == 0) write(*,*) (100.*c/(1.*grid % n_cells)), '% c_omplete...'  
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
  end do
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Interpolate density and velocity
    Us = grid % f(s) * U % n(c1) + (1.0 - grid % f(s)) * U % n(c2)
    Vs = grid % f(s) * V % n(c1) + (1.0 - grid % f(s)) * V % n(c2)
    Ws = grid % f(s) * W % n(c1) + (1.0 - grid % f(s)) * W % n(c2)
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
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
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

