!==============================================================================!
  subroutine ReaCom(restar) 
!------------------------------------------------------------------------------!
!   Reads second part of T-FlowS.cmn file.                                     ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical   :: restar
!----------------------------------[Calling]-----------------------------------!
  real      :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, j, l, m
  real              :: MresT, dummy
  character(len=80) :: answer, nammon
  real, allocatable :: mres(:), xm(:), ym(:), zm(:)
!==============================================================================!

  call Wait   

  ! The number of time steps
  if(this_proc  < 2) then 
    write(*,*) '# Enter the number of time steps: (',Ndt,') '
    write(*,*) '# (type 0 if you just want to analyse results)'
  end if
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp,*)  Ndt

  ! Starting time step for statistics 
  if(this_proc  < 2)  &
    write(*,*) '# Starting time step for statistics (',Nstat,') '
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp,*)  Nstat
  if(BUDG == YES) then
    read(inp(ts(2):te(2)),*) Nbudg
  end if

  if(this_proc  < 2)  & 
    write(*,*) '# Number of monitoring points:'
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp,*) Nmon 

  allocate (mres(Nmon))
  allocate (xm(Nmon))
  allocate (ym(Nmon))
  allocate (zm(Nmon))

  if(this_proc  < 2)  &
    write(*,*) '# Enter the coordinates of monitoring point(s)'
  do i=1,Nmon
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)  xm(i), ym(i), zm(i)
  end do

  ! Find the monitoring cells
  nammon=name 
  nammon(len_trim(name)+1:len_trim(name)+10)="-monit.000"
  l=len_trim(nammon) 
  do j=1,Nmon
    Mres(j)=HUGE
    do i=1,NC
      if(Distance(xm(j),ym(j),zm(j),  &
              xc(i),yc(i),zc(i))  < Mres(j)) then
        Cm(j)=i
        Mres(j)=Distance(xm(j),ym(j),zm(j),xc(i),yc(i),zc(i))
      end if
    end do
    MresT=Mres(j)
    call GloMin(MresT)
    if(MresT /= Mres(j)) then ! there is a cell which is nearer
      Cm(j) = 0               ! so erase this_proc monitoring point     
    end if 
  end do

  do j=1,Nmon
    if(Cm(j)  > 0) then
      if(j  <  10) then
        write(nammon(l  :l),'(I1)') j
      else if(j  < 100) then
        write(nammon(l-1:l),'(I2)') j
      else
        write(nammon(l-2:l),'(I3)') j
      end if
      if(Ndtt == 0) then
        open(10+j,file=nammon)
      else
        open(10+j,file=nammon,POSITION='APPEND')
      endif

      write(10+j,'(A24,3F16.6)')  &
            '# Monitoring point:',xc(Cm(j)),yc(Cm(j)),zc(Cm(j))
    end if
  end do

  ! Plane for calcution of overall mass fluxes
  do m=1,Nmat
    if(this_proc  < 2)  then
      write(*,*) '# Enter the coordinates of monitoring plane: (', &
                  xp(m), yp(m), zp(m), ' )'
    end if
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*) xp(m), yp(m), zp(m)
  end do

  ! Kind of simulation
  if(this_proc  < 2) then
    write(*,*) '# Type of simulation: '
    write(*,*) '# DNS      -> Direct Numerical Simulation'
    write(*,*) '# LES      -> Large Eddy Simulation'
    write(*,*) '# K_EPS    -> High Reynolds k-eps model.' 
    write(*,*) '# K_EPS_VV -> Durbin`s model.' 
    write(*,*) '# SPA_ALL  -> Spalart-Allmaras model.' 
    write(*,*) '# ZETA  -> k-eps-zeta-f model.' 
  endif
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'DNS') then
    SIMULA = DNS
  else if(answer == 'LES') then
    SIMULA = LES
  else if(answer == 'K_EPS') then
    SIMULA = K_EPS
  else if(answer == 'K_EPS_VV') then
    SIMULA = K_EPS_VV
  else if(answer == 'SPA_ALL') then
    SIMULA = SPA_ALL
  else if(answer == 'DES_SPA') then
    SIMULA = DES_SPA
  else if(answer == 'ZETA') then
    SIMULA = ZETA
  else if(answer == 'HYB_PITM') then
    SIMULA = HYB_PITM
  else if(answer == 'HYB_ZETA') then
    SIMULA = HYB_ZETA
  else if(answer == 'EBM') then
    SIMULA = EBM
  else if(answer == 'HJ') then
    SIMULA = HJ
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(SIMULA==EBM.or.SIMULA==HJ) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call To_Upper_Case(answer)
    if(answer == 'HYB') then
      MODE = HYB 
    else
      HYB = 10000
    end if
  end if 

  if(SIMULA==K_EPS) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call To_Upper_Case(answer)
    if(answer == 'HRE') then
      MODE = HRE
    else if(answer == 'LRE') then
      MODE = LRE 
    else
      if(this_proc  < 2) then
        write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                               cmn_line_count, ' Got a: ', answer
      endif
      stop
    end if
  end if

  if(SIMULA==LES.or.SIMULA==HYB_ZETA) then
    read(inp(ts(2):te(2)),'(A8)') answer
    call To_Upper_Case(answer)
    if(answer == 'SMAG') then
      MODE = SMAG 
    else if(answer == 'DYN') then
      MODE = DYN 
    else if(answer == 'WALE') then
      MODE = WALE
    else
      if(this_proc  < 2) then
        write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                               cmn_line_count, ' Got a: ', answer
      endif
      stop
    end if
  end if

  if(SIMULA  ==  LES.and.MODE == SMAG) then
    if(this_proc  < 2)  &
      write(*,*) '# C Smagorinsky = ', Cs0, ' enter the new value: '
    read(inp(ts(3):te(3)),*) Cs0
  endif 

    
  do m=1,Nmat
    if(SIMULA  ==  LES .or. SIMULA == DNS .or. SIMULA == DES_SPA &
      .or.SIMULA  ==  HYB_PITM .or. SIMULA == HYB_ZETA) then
      if(this_proc  < 2) then
        write(*,*) '# Do you want to shake the velocity field ?'
        write(*,*) '# YES -> shake'
        write(*,*) '# NO  -> don''t shake'
      end if
      call ReadC(CMN_FILE,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call To_Upper_Case(answer)
      if(answer == 'YES') then
        SHAKE(m) = YES
      else if(answer == 'NO') then
        SHAKE(m) = NO
      else
        if(this_proc  < 2) then
          write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                                 cmn_line_count, ' Got a: ', answer
        endif
        stop
      endif
      if(SHAKE(m) == YES) then
        if(this_proc < 2) &
          write(*,*) '# For how many time steps you want to shake ?'
        call ReadC(CMN_FILE,inp,tn,ts,te)
        read(inp,*) SHAKE_PER(m) 
        if(this_proc < 2) &
          write(*,*) '# Interval for shaking:', SHAKE_PER(m)
        call ReadC(CMN_FILE,inp,tn,ts,te)
        read(inp,*) SHAKE_INT(m)
      end if
    endif 
  end do

  if(.not. restar) call UnkAloc
  if(SIMULA  ==  K_EPS.or.SIMULA  ==  HYB_PITM) then
    Ce1 = 1.44
    Ce2 = 1.92
    Cmu = 0.09
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3 
    kappa = 0.41
    Elog  = 8.342
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
!    if(MODE == LRe) then
!      Ce1 = 1.55
!      Ce2 = 2.0
!    end if
  endif
 
  if(SIMULA  ==  K_EPS_VV) then
    Ce1   = 1.4
    Ce2   = 1.9
    Cmu   = 0.09
    CmuD  = 0.22
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    kappa = 0.41
    Elog  = 8.342
    Cl    = 0.22
    Ct    = 6.0
    Cni   = 85.0
    alpha = 0.045
    cf1   = 0.697
    cf2   = 4.4  
    cf3   = 1.6 
    Cf_1  = 1.4
    Cf_2  = 0.3
    Lim   = 11.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
    v_2 % Sigma = 1.0
  endif

  if(SIMULA  ==  EBM) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.15
    uu % Sigma = 1.0
    vv % Sigma = 1.0
    ww % Sigma = 1.0
    uv % Sigma = 1.0
    uw % Sigma = 1.0
    vw % Sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(SIMULA  ==  HJ) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.0 !1.15
    uu % Sigma = 1.0
    vv % Sigma = 1.0
    ww % Sigma = 1.0
    uv % Sigma = 1.0
    uw % Sigma = 1.0
    vw % Sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(SIMULA  ==  ZETA.or.SIMULA  ==  HYB_ZETA) then
    Ce1   = 1.4
    Ce2   = 1.9
    Cmu   = 0.09 
    CmuD  = 0.22
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    kappa = 0.41
    Elog  = 8.342
    Cl    = 0.36
    Ct    = 6.0
    Cni   = 85.0
    alpha = 0.012
    cf1   = 0.697
    cf2   = 4.4
    cf3   = 1.6
    Cf_1  = 1.4
    Cf_2  = 0.3
    Lim   = 11.0
    Kin % Sigma = 1.0
    Eps % Sigma = 1.3
    v_2 % Sigma = 1.2
  endif

  if(SIMULA == SPA_ALL .or. SIMULA==DES_SPA) then
    kappa  = 0.41
    Cb1    = 0.1355
    Cb2    = 0.622
    Cvis1  = 7.1
    Cw2    = 0.3
    Cw3    = 2.0
    VIS % Sigma = 2.0/3.0
    Cw1    = Cb1/kappa**2.0 + (1+Cb2)/VIS % Sigma
    SIGMAv = 2.0/3.0
  end if

  ! Time stepping scheme
  if(this_proc  < 2) then
    write(*,*) '# Algorythm for time-integration: '
    write(*,*) '# SIMPLE [Nini] -> S. I. M. P. L. E.'
    write(*,*) '# FRACTION      -> Fractional step method'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'SIMPLE') then
    ALGOR = SIMPLE
    Nini  = 10
    if(tn==2) read(inp(ts(2):te(2)),*) Nini
  else if(answer == 'FRACTION') then
    ALGOR = FRACT
    Nini  = 1
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(ALGOR == SIMPLE) then
    if(this_proc < 2) write(*,*) '# Under Relaxation Factor for velocity (',U % URF,')'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)  U % URF
    if(this_proc < 2) write(*,*) '# Under Relaxation Factor for pressure (',P % URF,')'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)  P % URF
    if(HOT == YES) then
      if(this_proc < 2) write(*,*) '# Under Relaxation Factor for temperature (',T % URF,')'
      call ReadC(CMN_FILE,inp,tn,ts,te)
      read(inp,*)  T % URF
    end if
    if(SIMULA /= LES .and. SIMULA /= DNS) then
      if(this_proc < 2) write(*,*) '# Under Relaxation Factor for turbulent variables (',T % URF,')'
      call ReadC(CMN_FILE,inp,tn,ts,te)
      read(inp,*)  URFT
    end if
  endif
  Kin % URF   = URFT
  Eps % URF   = URFT
  v_2 % URF   = URFT
  f22 % URF   = URFT 
  VIS % URF   = URFT
  uu  % URF   = URFT
  vv  % URF   = URFT
  ww  % URF   = URFT
  uv  % URF   = URFT
  uw  % URF   = URFT
  vw  % URF   = URFT


  if(this_proc  < 2) then
    write(*,*) '# Integration of inertial terms: '
    write(*,*) '# LIN -> Linear'
    write(*,*) '# PAR -> Parabolic'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'LIN') then
    INERT = LIN
  else if(answer == 'PAR') then
    INERT = PAR
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(this_proc  < 2) then
    write(*,*) '# Integration of convective terms: '
    write(*,*) '# AB -> Adams-Bashforth'
    write(*,*) '# CN -> Crank-Nicholson'
    write(*,*) '# FI -> Fully Implicit'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'AB') then
    CONVEC = AB
  else if(answer == 'CN') then
    CONVEC = CN
  else if(answer == 'FI') then
    CONVEC = FI
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(this_proc  < 2) then
    write(*,*) '# Integration of diffusive terms: '
    write(*,*) '# AB -> Adams-Bashforth'
    write(*,*) '# CN -> Crank-Nicholson'
    write(*,*) '# FI -> Fully Implicit'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'AB') then
    DIFFUS = AB
  else if(answer == 'CN') then
    DIFFUS = CN
  else if(answer == 'FI') then
    DIFFUS = FI
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(this_proc  < 2) then
    write(*,*) '# Integration of cross-diffusive terms: '
    write(*,*) '# AB -> Adams-Bashforth'
    write(*,*) '# CN -> Crank-Nicholson'
    write(*,*) '# FI -> Fully Implicit'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'AB') then
    CROSS = AB
  else if(answer == 'CN') then
    CROSS = CN
  else if(answer == 'FI') then
    CROSS = FI
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  ! Upwind blending
  do m=1,Nmat
    URFC(m) = 1.0
    if(this_proc  < 2) then
      write(*,*) '# Convetive schemes for momentum equation:'
      write(*,*) '# Do you want to use upwind blending: '
      write(*,*) '# YES       -> use blening'
      write(*,*) '# NO        -> don''t use blending'
      write(*,*) '# CDS       -> central differencing'
      write(*,*) '# LUDS      -> linear upwind'
      write(*,*) '# QUICK     -> self descriptive'
      write(*,*) '# MINMOD    -> self descriptive'
      write(*,*) '# SMART     -> self descriptive'
      write(*,*) '# AVL_SMART -> self descriptive'
    endif 
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp(ts(1):te(1)),'(A)')  answer
    call To_Upper_Case(answer)
    if(answer == 'BLEND_CDS_UDS') then
      BLEND(m) = YES
      if(tn==2) read(inp(ts(2):te(2)),*) URFC(m)
    else if(answer == 'NO') then
      BLEND(m) = NO
    else if(answer == 'UDS') then
      BLEND(m) = YES
      URFC(m)  = 0.0
    else if(answer == 'CDS') then
      BLEND(m) = CDS 
    else if(answer == 'LUDS') then
      BLEND(m) = LUDS
    else if(answer == 'QUICK') then
      BLEND(m) = QUICK
    else if(answer == 'MINMOD') then
      BLEND(m) = MINMOD
    else if(answer == 'SMART') then
      BLEND(m) = SMART 
    else if(answer == 'AVL_SMART') then
      BLEND(m) = AVL_SMART 
    else if(answer == 'SUPERBEE') then
      BLEND(m) = SUPERBEE 
    else if(answer == 'GAMMA') then
      BLEND(m) = GAMMA 
    else
      if(this_proc  < 2) then
        write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                               cmn_line_count, ' Got a: ', answer
      endif
      stop
    endif
  end do

  if(HOT==YES) then
    do m=1,Nmat
      URFC_Tem(m) = 1.0
      if(this_proc  < 2) then
        write(*,*) '# Convetive schemes for energy equation:'
      endif 
      call ReadC(CMN_FILE,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call To_Upper_Case(answer)
      if(answer == 'BLEND_TEM_CDS_UDS') then
        BLEND_TEM(m) = YES
        if(tn==2) read(inp(ts(2):te(2)),*) URFC_Tem(m)
      else if(answer == 'NO') then
        BLEND_TEM(m) = NO
      else if(answer == 'UDS') then
        BLEND_TEM(m) = YES
        URFC_Tem(m)  = 0.5
      else if(answer == 'CDS') then
        BLEND_TEM(m) = CDS 
      else if(answer == 'LUDS') then
        BLEND_TEM(m) = LUDS
      else if(answer == 'QUICK') then
        BLEND_TEM(m) = QUICK
      else if(answer == 'MINMOD') then
        BLEND_TEM(m) = MINMOD
      else if(answer == 'SMART') then
        BLEND_TEM(m) = SMART 
      else if(answer == 'AVL_SMART') then
        BLEND_TEM(m) = AVL_SMART 
      else if(answer == 'SUPERBEE') then
        BLEND_TEM(m) = SUPERBEE 
      else if(answer == 'GAMMA') then
        BLEND_TEM(m) = GAMMA 
      else
        if(this_proc  < 2) then
          write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                                 cmn_line_count, ' Got a: ', answer
        endif
        stop
      endif
    end do
  end if

  if(SIMULA/=LES.and.SIMULA/=DNS) then
    do m=1,Nmat
      URFC_Tur(m) = 1.0
      if(this_proc  < 2) then
        write(*,*) '# Convetive schemes for transport equation:'
      endif 
      call ReadC(CMN_FILE,inp,tn,ts,te)
      read(inp(ts(1):te(1)),'(A)')  answer
      call To_Upper_Case(answer)
      if(answer == 'BLEND_TUR_CDS_UDS') then
        BLEND_TUR(m) = YES
        if(tn==2) read(inp(ts(2):te(2)),*) URFC_Tur(m)
      else if(answer == 'NO') then
        BLEND_TUR(m) = NO
      else if(answer == 'UDS') then
        BLEND_TUR(m) = YES
        URFC_Tur(m)  = 0.0
      else if(answer == 'CDS') then
        BLEND_TUR(m) = CDS 
      else if(answer == 'LUDS') then
        BLEND_TUR(m) = LUDS
      else if(answer == 'QUICK') then
        BLEND_TUR(m) = QUICK
      else if(answer == 'MINMOD') then
        BLEND_TUR(m) = MINMOD
      else if(answer == 'SMART') then
        BLEND_TUR(m) = SMART 
      else if(answer == 'AVL_SMART') then
        BLEND_TUR(m) = AVL_SMART 
      else if(answer == 'SUPERBEE') then
        BLEND_TUR(m) = SUPERBEE 
      else if(answer == 'GAMMA') then
        BLEND(m) = GAMMA 
      else
        if(this_proc  < 2) then
          write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                                 cmn_line_count, ' Got a: ', answer
        endif
        stop
      endif
    end do
  end if

  ! Solver parameters
  if(this_proc  < 2) then
    write(*,*) '# Preconditioning of the system matrix: '
    write(*,*) '# NO -> No preconditioning'
    write(*,*) '# DI -> Diagonal preconditioning'
    write(*,*) '# IC -> Incomplete Cholesky'
  endif 
  call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp(ts(1):te(1)),'(A)')  answer
  call To_Upper_Case(answer)
  if(answer == 'NO') then
    PREC = 0 
  else if(answer == 'DI') then
    PREC = 1 
  else if(answer == 'IC') then
    PREC = 2 
  else
    if(this_proc  < 2) then
      write(*,'(A,I3,A,A)') 'Error in T-FlowS.cmn file in line ', &
                             cmn_line_count, ' Got a: ', answer
    endif
    stop
  endif

  if(this_proc  < 2)  &
    write(*,*) '# Tolerance for velocity solver: (',U % STol,' )'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)    U % STol
    V % Stol     = U % Stol
    W % Stol     = U % Stol
  if(this_proc  < 2)  &
    write(*,*) '# Tolerance for pressure solver: (',PP % STol,' )'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)   PP % STol
    P % Stol = PP % Stol
  if(SIMULA/=LES.and.SIMULA/=DNS) then
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)    Kin % STol
    Eps % Stol   = Kin % Stol
    v_2 % Stol   = Kin % Stol
    f22 % Stol   = Kin % Stol
    VIS % Stol   = Kin % Stol
    uu  % Stol   = Kin % Stol
    vv  % Stol   = Kin % Stol
    ww  % Stol   = Kin % Stol
    uv  % Stol   = Kin % Stol
    uw  % Stol   = Kin % Stol
    vw  % Stol   = Kin % Stol
  end if
  if(HOT == YES) then
    if(this_proc  < 2)  &
      write(*,*) '# Tolerance for temperature solver: (',T % STol,' )'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)    T % STol
  end if
 
  if(ALGOR == SIMPLE) then
    if(this_proc  < 2)  &
      write(*,*) '# Tolerance for SIMPLE: (',SIMTol,' )'
    call ReadC(CMN_FILE,inp,tn,ts,te)
    read(inp,*)   SIMTol
  endif     

  ! Time step
  if(this_proc  < 2)  &
    write(*,*) '# Time step: (',dt,' )'
    call ReadC(CMN_FILE,inp,tn,ts,te)
  read(inp,*)   dt

  ! Wall velocity 
  do m=1,Nmat
    if(this_proc  < 2)  &
      write(*,*) '# Enter Pdrop (x, y, z) for domain ', m
    call ReadC(CMN_FILE,inp,tn,ts,te) 
    if(.not. restar) then 
      read(inp,*)  PdropX(m), PdropY(m), PdropZ(m)
      UTau(m) = sqrt(abs(PdropX(m))) ! delta=1, nu=1 
      VTau(m) = sqrt(abs(PdropY(m))) ! delta=1, nu=1 
      WTau(m) = sqrt(abs(PdropZ(m))) ! delta=1, nu=1 
    else
      read(inp,*)  dummy, dummy, dummy 
    end if     
  end do

  ! Mass fluxes
  do m=1,Nmat
    if(this_proc  < 2) then
      write(*,*) '# Enter the wanted mass flux through domain ', m
      write(*,*) '# (type 0.0 to keep the pressure drop constant)'
    endif 
    call ReadC(CMN_FILE,inp,tn,ts,te) 
    if(.not. restar) read(inp,*)  FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    if(restar)       read(inp,*)  dummy, dummy, dummy 
  end do

  call Wait   

  end subroutine ReaCom

