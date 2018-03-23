!==============================================================================!
  subroutine ReaCom(grid, restar)
!------------------------------------------------------------------------------!
!   Reads second part of T-FlowS.cmn file.                                     ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use les_mod
  use Comm_Mod
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  logical           :: restar
!----------------------------------[Calling]-----------------------------------!
  real              :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, j, l, m
  real              :: MresT
  character(len=80) :: nammon
  real, allocatable :: mres(:), xm(:), ym(:), zm(:)
!==============================================================================!

  call Comm_Mod_Wait

  ! Retreive monitoring points
  call Control_Mod_Monitoring_Points(Nmon, xm, ym, zm)

  allocate(Mres(Nmon))

  ! Find the monitoring cells
  nammon=problem_name
  nammon(len_trim(problem_name)+1:len_trim(problem_name)+10)='-monit.000'
  l=len_trim(nammon)
  do j = 1, Nmon
    Mres(j)=HUGE
    do i = 1, grid % n_cells
      if(Distance(xm(j),ym(j),zm(j),  &
              grid % xc(i),grid % yc(i),grid % zc(i))  < Mres(j)) then
        Cm(j)=i
        Mres(j)=Distance(xm(j),ym(j),zm(j),grid % xc(i),grid % yc(i),grid % zc(i))
      end if
    end do
    MresT=Mres(j)
    call Comm_Mod_Global_Min_Real(MresT)
    if(MresT /= Mres(j)) then ! there is a cell which is nearer
      Cm(j) = 0               ! so erase this_proc monitoring point
    end if
  end do

  do j = 1, Nmon
    if(Cm(j)  > 0) then
      if(j  <  10) then
        write(nammon(l  :l),'(I1)') j
      else if(j  < 100) then
        write(nammon(l-1:l),'(I2)') j
      else
        write(nammon(l-2:l),'(I3)') j
      end if
      if(.not. restar) then
        open(10+j,file=nammon)
      else
        open(10+j,file=nammon,POSITION='APPEND')
      endif

      write(10+j,'(A24,3F16.6)')  &
            '# Monitoring point:',grid % xc(Cm(j)),grid % yc(Cm(j)),grid % zc(Cm(j))
    end if
  end do

  ! Plane for calcution of overall mass fluxes
  do m = 1, grid % n_materials
    call Control_Mod_Point_For_Monitoring_Plane(bulk(m) % xp,  &
                                                bulk(m) % yp,  &
                                                bulk(m) % zp)
  end do

  if(this_proc < 2) &
    print *, 'Calling for gravitationa vector'
  call Control_Mod_Gravitational_Vector(grav_x, grav_y, grav_z)

  if(.not. restar) call Allocate_Variables(grid)

  if(this_proc < 2) &
    print *, 'Calling for turbulence model'
  call Control_Mod_Turbulence_Model(.true.)

  if(turbulence_model == K_EPS .or.  &
     turbulence_model == HYBRID_PITM) then
    Ce1 = 1.44
    Ce2 = 1.92
    Cmu = 0.09
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    kappa = 0.41
    Elog  = 8.342
    Kin % sigma = 1.0
    Eps % sigma = 1.3
!    if(MODE == LRe) then
!      Ce1 = 1.55
!      Ce2 = 2.0
!    end if
  endif

  if(turbulence_model == REYNOLDS_STRESS_MODEL) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    kin % sigma = 1.0
    eps % sigma = 1.15
    uu % sigma = 1.0
    vv % sigma = 1.0
    ww % sigma = 1.0
    uv % sigma = 1.0
    uw % sigma = 1.0
    vw % sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(turbulence_model == HANJALIC_JAKIRLIC) then
    Ce1   = 1.44
    Ce2   = 1.83
    Ce3   = 0.55
    CmuD  = 0.21
    Cmu25 = sqrt(sqrt(Cmu))
    Cmu75 = Cmu25**3
    Cl    = 0.161
    Ct    = 6.0
    Cni   = 80.0
    kin % sigma = 1.0
    eps % sigma = 1.0 !1.15
    uu % sigma = 1.0
    vv % sigma = 1.0
    ww % sigma = 1.0
    uv % sigma = 1.0
    uw % sigma = 1.0
    vw % sigma = 1.0
    g1    = 3.4
    g1_star = 1.8
    g2    = 4.2
    g3    = 0.8
    g3_star = 1.3
    g4    = 1.25
    g5    = 0.4
  endif

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
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
    kin  % sigma = 1.0
    eps  % sigma = 1.3
    zeta % sigma = 1.2
  endif

  if(turbulence_model == SPALART_ALLMARAS .or.  &
     turbulence_model == DES_SPALART) then
    kappa  = 0.41
    Cb1    = 0.1355
    Cb2    = 0.622
    Cvis1  = 7.1
    Cw2    = 0.3
    Cw3    = 2.0
    VIS % sigma = 2.0/3.0
    Cw1    = Cb1/kappa**2.0 + (1+Cb2)/VIS % sigma
    SIGMAv = 2.0/3.0
  end if

  ! Wall velocity
  if(.not. restar) then
    do m=1,grid % n_materials
      call Control_Mod_Mass_Flow_Rates(bulk(m) % p_drop_x,  &
                                       bulk(m) % p_drop_y,  &
                                       bulk(m) % p_drop_z)
    end do
  end if

  ! Mass fluxes
  if(.not. restar) then
    do m=1,grid % n_materials
      call Control_Mod_Mass_Flow_Rates(bulk(m) % flux_x_o,  &
                                       bulk(m) % flux_y_o,  &
                                       bulk(m) % flux_z_o)
    end do
  end if

  call Comm_Mod_Wait

  end subroutine