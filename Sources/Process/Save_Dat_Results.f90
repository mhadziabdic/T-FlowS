!==============================================================================!
  subroutine Save_Dat_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   Writes: name.dat                                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use par_mod, only: this_proc
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Constants_Pro_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c,  c2,  n, s, Nadd
  integer           :: Nfac(10), NtotFac
  character(len=80) :: name_out, store_name
!==============================================================================!

  ! Store the name
  store_name = problem_name     

  problem_name = name_save  

  call Wait 

  !----------------------!
  !                      !
  !   Create .dat file   !
  !                      !
  !----------------------!
  call Name_File(this_proc, name_out, '.dat', len_trim('.dat'))
  open(9, file=name_out)
  if(this_proc  < 2) print *, '# Now creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(A34)') '(0 "============================")'
  write(9,'(A34)') '(0 "Created by T-Rex - Processor")'
  write(9,'(A34)') '(0 "============================")'
  Nadd = 0
  if(SIMULA == LES .or. SIMULA == DES_SPA) then
    Nadd = Nadd + 9   ! Um, Vm, Wm, uu, vv, ww, uv, uw, vw
  end if
  if(HOT == YES) then
    Nadd = Nadd + 1
    if(SIMULA == LES .or. SIMULA == DES_SPA) then
      Nadd = Nadd + 5 ! Tm, TT, uT, vT, wT
    end if
  endif
  write(9,'(A3, I3, A2)') '(0 ', Nadd, ' )'

  !-------------!
  !   Results   !
  !-------------+---------------------------------------------------!
  !   (300 (sub_section_id  zone_id  size  n_time_levels  n_phases  !
  !         first_id  last_id)                                      !         
  !-----------------------------------------------------------------!

  !----------------!
  !   Velocities   !
  !----------------!
  write(9,'(A16)') '(0 "velocities")'
  write(9,'(A16,I9,I9,A2)') '(300 (2 1 3 0 0 ',  1, grid % n_cells, ')(' 
  do c = 1, grid % n_cells
    write(9,'(3F14.6)') U % n(c), V % n(c), W % n(c) 
  end do  
  write(9,'(A2)') '))'

  !--------------!
  !   Pressure   !
  !--------------!
  call Save_Dat_Scalar(grid, 'pressure', 1, P % n)

  !-----------------------!
  !   Aditional scalars   !
  !-----------------------!
  Nadd = 0 

  !-------------------!
  !   Other scalars   !
  !-------------------!
  if(HOT == YES) then
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is temperature' 
    call Save_Dat_Scalar('temperature', 699+Nadd, T % n)
  end if
  
  if(SIMULA == K_EPS) then
    ! VIS  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is K' 
    call Save_Dat_Scalar('Kin', 699+Nadd, Kin % n)
    ! VISt  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Eps' 
    call Save_Dat_Scalar('Eps', 699+Nadd, Eps % n)
  end if

  if(SIMULA == K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
    ! VIS  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is K' 
    call Save_Dat_Scalar('Kin', 699+Nadd, Kin % n)
    ! VISt  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Eps' 
    call Save_Dat_Scalar('Eps', 699+Nadd, Eps % n)
    ! Vort  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is v_2' 
    call Save_Dat_Scalar('v_2', 699+Nadd, v_2 % n)
    ! WallDs
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is f22' 
    call Save_Dat_Scalar('f22', 699+Nadd, f22 % n)

  end if ! SIMULA=K_EPS_VV

  if(SIMULA == SPA_ALL) then
    ! VIS  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is VIS' 
    call Save_Dat_Scalar('VIS', 699+Nadd, VIS % n)
    ! VISt  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is VISt' 
    call Save_Dat_Scalar('VISt', 699+Nadd, VISt)
    ! Vort  
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Vort' 
    call Save_Dat_Scalar('Vort', 699+Nadd, Vort)
    ! WallDs
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is WallDs' 
    call Save_Dat_Scalar('WallDs', 699+Nadd, WallDs)
  end if ! SIMULA=SPA_ALL

  if(SIMULA == LES .or. SIMULA == DES_SPA) then
    ! Umean
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Umean' 
    call Save_Dat_Scalar('Umean', 699+Nadd, U % mean)
    ! Vmean
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Vmean' 
    call Save_Dat_Scalar('Vmean', 699+Nadd, V % mean)
    ! Wmean
    Nadd = Nadd + 1
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Wmean' 
    call Save_Dat_Scalar('Wmean', 699+Nadd, W % mean)
    ! uu
    Nadd = Nadd + 1
    PP % n = uu % mean - U % mean * U % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uu' 
    call Save_Dat_Scalar('uu', 699+Nadd, PP % n)
    ! vv
    Nadd = Nadd + 1
    PP % n = vv % mean - V % mean * V % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is vv' 
    call Save_Dat_Scalar('vv', 699+Nadd, PP % n) 
    ! ww
    Nadd = Nadd + 1
    PP % n = ww % mean - W % mean * W % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is ww' 
    call Save_Dat_Scalar('ww', 699+Nadd, PP % n)
    ! uv
    Nadd = Nadd + 1
    PP % n = uv % mean - U % mean * V % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uv' 
    call Save_Dat_Scalar('uv', 699+Nadd, PP % n)
    ! uw
    Nadd = Nadd + 1
    PP % n = uw % mean - U % mean * W % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uw' 
    call Save_Dat_Scalar('uw', 699+Nadd, PP % n) 
    ! vw
    Nadd = Nadd + 1
    PP % n = vw % mean - V % mean * W % mean
    if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is ww' 
    call Save_Dat_Scalar('vw', 699+Nadd, PP % n)

    if(HOT == YES) then
      ! Tmean
      Nadd = Nadd + 1
      if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is Tmean' 
      call Save_Dat_Scalar('Tmean', 699+Nadd, T % mean)
      ! TT    
      Nadd = Nadd + 1
      PP % n = TT % mean - T % mean * T % mean
      if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is TT' 
      call Save_Dat_Scalar('TT', 699+Nadd, PP % n)
      ! uT    
      Nadd = Nadd + 1
      PP % n = uT % mean - U % mean * T % mean
      if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is uT' 
      call Save_Dat_Scalar('uT', 699+Nadd, PP % n)
      ! vT    
      Nadd = Nadd + 1
      PP % n = vT % mean - V % mean * T % mean
      if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is vT' 
      call Save_Dat_Scalar('vT', 699+Nadd, PP % n)
      ! wT    
      Nadd = Nadd + 1
      PP % n = wT % mean - W % mean * T % mean
      if(this_proc  < 2) write(*,'(A,I3,A)') '# Scalar ', Nadd, ' is wT' 
      call Save_Dat_Scalar('wT', 699+Nadd, PP % n)
    end if
  end if ! SIMULA == LES .or. SIMULA == DES_SPA

  !-----------------------------------!
  !   This is important for connect   !
  !-----------------------------------!
  NtotFac = 0
  do n=1,10   ! Browse through boundary condition types
    Nfac(n) = 0
    do s = 1, grid % n_faces   ! Count the faces with boundary condition "n"
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == n) Nfac(n)=Nfac(n)+1
      end if
    end do    ! sides

    ! Prepare for next boundary
    if(this_proc  < 2) print *, 'Number of faces:', Nfac(n), NtotFac+1, NtotFac+Nfac(n)
    NtotFac = NtotFac+Nfac(n)
  end do   ! n -> boundary condition types

  write(9,'(I8)') grid % n_cells
  do n=1,10
    write(9,'(I8)') Nfac(n)
  end do

  close(9)

  ! Restore the name
  problem_name = store_name  

  end subroutine
