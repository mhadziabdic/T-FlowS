!==============================================================================!
  subroutine ReaCom(grid, restar)
!------------------------------------------------------------------------------!
!   This subroutine used to read the second part of T-FlowS.cmn file.          ! 
!   Now, with T-FlowS.cmn file gone, it is pretty much oboslete and candidate  !
!   for deletion.  At the moment, it is merely defining monitoring points and  !
!   initializes constants in turbulence models.                                !
!   The latter is clearly miss-placed here.                                    !
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
    if(abs(MresT - Mres(j)) <= TINY) then ! there is a cell which is nearer
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

  if(turbulence_model == K_EPS) then
    call Constants_K_Eps()
  endif

  if(turbulence_model == REYNOLDS_STRESS_MODEL) then
    call Constants_Reynolds_Stress()
  endif

  if(turbulence_model == HANJALIC_JAKIRLIC) then
    call Constants_Hanjalic_Jakirlic()
  endif

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    call Constants_K_Eps_Zeta_F()
  endif

  if(turbulence_model == SPALART_ALLMARAS .or.  &
     turbulence_model == DES_SPALART) then
    call Constants_Spalart_Allmaras()
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
