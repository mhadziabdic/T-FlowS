!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured Finite Volume 'LES'/RANS solver.                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use les_mod
  use Comm_Mod
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Grad_Mod
  use Bulk_Mod
  use Var_Mod
  use Solvers_Mod, only: D
  use Info_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Calling]------------------------------------!
  real :: Correct_Velocity
!----------------------------------[Locals]------------------------------------!
  integer           :: i, m, n
  real              :: mass_res, wall_time_start, wall_time_current
  character(len=80) :: name_save
  logical           :: restar, multiple, save_now, exit_now
  type(Grid_Type)   :: grid        ! grid used in computations
  real              :: time        ! physical time
  real              :: dt          ! time step
  integer           :: first_dt    ! first time step in this run
  integer           :: last_dt     ! number of time steps
  integer           :: max_ini     ! max number of inner iterations
  integer           :: min_ini     ! min number of inner iterations
  integer           :: n_stat      ! starting time step for statistic
  integer           :: ini         ! inner iteration counter
  integer           :: bsi, rsi    ! backup and results save interval
  real              :: simple_tol  ! tolerance for SIMPLE algorithm
  character(len=80) :: coupling    ! pressure velocity coupling
!--------------------------------[Interfaces]----------------------------------!
  interface
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine UserProbe2D(namAut)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      use Name_Mod
      use Flow_Mod
      implicit none
      character, optional :: namAut*(*)
    end subroutine
  end interface
!==============================================================================!

  ! Get starting time
  call cpu_time(wall_time_start)

  !------------------------------!
  !   Start parallel execution   !
  !------------------------------!
  call Comm_Mod_Start

  !--------------------------------!
  !   Splash out the logo screen   !
  !--------------------------------!
  if(this_proc  < 2) then
    call Logo_Pro
  endif

  !---------------------------------------------!
  !   Open control file and read problem name   !
  !---------------------------------------------!
  call Control_Mod_Open_File()
  call Control_Mod_Problem_Name(problem_name)

  ! Initialize parameters -> bad practice, candidate for deletion
  call IniPar

  ! Load the finite volume grid
  call Load_Cns(grid, this_proc)

  ! Read problem type -> bad practice, candidate for deletion
  ! call Read_Problem

  call Allocate_Memory(grid)
  call Load_Geo       (grid, this_proc)
  call Comm_Mod_Load_Buffers
  call Comm_Mod_Exchange(grid, grid % vol(-grid % n_bnd_cells))

  call Matrix_Mod_Topology(grid, A)
  call Matrix_Mod_Topology(grid, D)

  call Comm_Mod_Wait

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!

  ! Get the number of time steps from the control file
  call Control_Mod_Number_Of_Time_Steps(last_dt, verbose=.true.)
  call Control_Mod_Starting_Time_Step_For_Statistics(n_stat, verbose=.true.)

  ! First time step is one, unless read from restart otherwise
  first_dt = 0
  call Load_Restart(grid, first_dt, restar)

  if(restar) then
    call Load_Boundary_Conditions(grid, .false.)
    call Load_Physical_Properties(grid)
  end if

  ! Read command file (T-FlowS.cmn)
  call ReaCom(grid, restar)

  ! Initialize variables
  if(.not. restar) then
    call Load_Boundary_Conditions(grid, .true.)
    call Load_Physical_Properties(grid)
    call Initialize_Variables(grid)
    call Comm_Mod_Wait
  end if

  ! Interpolate between diff. meshes
  call Load_Ini(grid)

  ! Check if there are more materials
  multiple = .false.
  i = StateMat(1)
  do m = 1, grid % n_materials
    if(StateMat(m) /= i) multiple = .true.
  end do

  ! Loading data from previous computation
  !  if(this_proc<2)  &
  !    print *,'Reading data from previous computation on the same mesh'
  call Load_Restart_Ini(grid)

  ! Prepare ...
  call Calculate_Face_Geometry(grid)
  call Bulk_Mod_Monitoring_Planes_Areas(grid, bulk)
  call Grad_Mod_Find_Bad_Cells         (grid)
  if(turbulence_model == LES                 .and.  &
     turbulence_model_variant == SMAGORINSKY .and.  &
     .not. restar)                                  &
     call Find_Nearest_Wall_Cell(grid)

  ! Prepare the gradient matrix for velocities
  call Compute_Gradient_Matrix(grid, .true.)

  ! Prepare matrix for fractional step method
  call Control_Mod_Pressure_Momentum_Coupling(coupling)
  if(coupling == 'PROJECTION') then
    call Pressure_Matrix_Fractional(grid, dt)
  end if

  ! Print the areas of monitoring planes
  if(this_proc < 2) then
    do m = 1, grid % n_materials
      write(*,'(a6,i2,a2,1pe12.3)') ' # Ax(',m,')=', bulk(m) % area_x
      write(*,'(a6,i2,a2,1pe12.3)') ' # Ay(',m,')=', bulk(m) % area_y
      write(*,'(a6,i2,a2,1pe12.3)') ' # Az(',m,')=', bulk(m) % area_z
    end do
  end if

  !---------------!
  !               !
  !   Time loop   !
  !               !
  !---------------!

  call Control_Mod_Time_Step(dt, verbose=.true.)
  call Control_Mod_Backup_Save_Interval(bsi, verbose=.true.)
  call Control_Mod_Results_Save_Interval(rsi, verbose=.true.)


  ! new feature: determine if code was compiles with CGNS
  COMPILED_WITH_CGNS = .true. ! assume it was
  call Save_Results(grid, "mesh")

  do n = first_dt + 1, last_dt

    time = time + dt

    ! Start info boxes.
    call Info_Mod_Time_Start()
    call Info_Mod_Iter_Start()
    call Info_Mod_Bulk_Start()

    ! Initialize and print time info box
    call cpu_time(wall_time_current)
    call Info_Mod_Time_Fill( n, time, (wall_time_current-wall_time_start) )
    call Info_Mod_Time_Print()

    if(turbulence_model == DES_SPALART) then
      call Calculate_Shear_And_Vorticity(grid)
      call Calculate_Vorticity (grid, u % n, v % n, w % n, vort)
    end if

    if(turbulence_model == LES) then
      call Calculate_Shear_And_Vorticity(grid)
      if(turbulence_model_variant == DYNAMIC) call Calculate_Sgs_Dynamic(grid)
      if(turbulence_model_variant == WALE)    call Calculate_Sgs_Wale(grid)
      call Calculate_Sgs(grid)
    end if

    If(turbulence_model == HYBRID_K_EPS_ZETA_F) then
      call Calculate_Sgs_Dynamic(grid)
      call Calculate_Sgs_Hybrid(grid)
    end if

    call Convective_Outflow(grid, dt)
    if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
       turbulence_model == HANJALIC_JAKIRLIC)       &
      call Calculate_Vis_T_Rsm(grid)

    !--------------------------!
    !   Inner-iteration loop   !
    !--------------------------!
    if(coupling == 'PROJECTION') then
      max_ini = 1
    else
      call Control_Mod_Max_Simple_Iterations(max_ini)
      call Control_Mod_Min_Simple_Iterations(min_ini)
    end if

    do ini = 1, max_ini  !  PROJECTION & SIMPLE

      call Info_Mod_Iter_Fill(ini)

      call Grad_Mod_For_P(grid, p % n, p % x, p % y, p % z)

      ! Compute velocity gradients
      call Grad_Mod_For_Phi(grid, u % n, 1, u % x, .true.)
      call Grad_Mod_For_Phi(grid, u % n, 2, u % y, .true.)
      call Grad_Mod_For_Phi(grid, u % n, 3, u % z, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 1, v % x, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 2, v % y, .true.)
      call Grad_Mod_For_Phi(grid, v % n, 3, v % z, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 1, w % x, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 2, w % y, .true.)
      call Grad_Mod_For_Phi(grid, w % n, 3, w % z, .true.)

      ! u velocity component
      call Compute_Momentum(grid, dt, ini, u,          &
                  u % x,   u % y,   u % z,             &
                  grid % sx,   grid % sy,   grid % sz, &
                  grid % dx,   grid % dy,   grid % dz, &
                  p % x,   v % x,   w % x)      ! dP/dx, dV/dx, dW/dx

      ! v velocity component
      call Compute_Momentum(grid, dt, ini, v,          &
                  v % y,   v % x,   v % z,             &
                  grid % sy,   grid % sx,   grid % sz, &
                  grid % dy,   grid % dx,   grid % dz, &
                  p % y,   u % y,   w % y)      ! dP/dy, dU/dy, dW/dy

      ! w velocity component
      call Compute_Momentum(grid, dt, ini, w,          &
                  w % z,   w % x,   w % y,             &
                  grid % sz,   grid % sx,   grid % sy, &
                  grid % dz,   grid % dx,   grid % dy, &
                  p % z,   u % z,   v % z)      ! dP/dz, dU/dz, dV/dz

      if(coupling == 'PROJECTION') then
        call Comm_Mod_Exchange(grid, A % sav)
        call Balance_Mass(grid)
        call Compute_Pressure_Fractional(grid, dt)
      endif
      if(coupling == 'SIMPLE') then
        call Comm_Mod_Exchange(grid, A % sav)
        call Balance_Mass(grid)
        call Compute_Pressure_Simple(grid)
      end if

      call Grad_Mod_For_P(grid,  pp % n, p % x, p % y, p % z)

      call Bulk_Mod_Compute_Fluxes(grid, bulk, flux)
      mass_res = Correct_Velocity(grid, dt) !  project the velocities

      ! Temperature
      if(heat_transfer == YES) then
        call Compute_Temperature(grid, dt, ini, t)
      end if

      ! Rans models
      if(turbulence_model == K_EPS) then 

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call Calculate_Shear_And_Vorticity(grid)

        call Compute_Turbulent(grid, dt, ini, kin, n)
        call Compute_Turbulent(grid, dt, ini, eps, n)

        call Calculate_Vis_T_K_Eps(grid)

        if(heat_transfer == YES) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model == K_EPS_ZETA_F     .or.  &
         turbulence_model == HYBRID_K_EPS_ZETA_F) then
        call Calculate_Shear_And_Vorticity(grid)

        call Compute_Turbulent(grid, dt, ini, kin, n)
        call Compute_Turbulent(grid, dt, ini, eps, n)
        call Update_Boundary_Values(grid)

        call Compute_F22(grid, ini, f22)
        call Compute_Turbulent(grid, dt, ini, zeta, n)

        call Calculate_Vis_T_K_Eps_Zeta_F(grid)

        if(heat_transfer == YES) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
         turbulence_model == HANJALIC_JAKIRLIC) then

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        if(turbulence_model == REYNOLDS_STRESS_MODEL) then
          call Time_And_Length_Scale(grid)
        end if

        call Grad_Mod_For_Phi(grid, u % n, 1, u % x,.true.)    ! dU/dx
        call Grad_Mod_For_Phi(grid, u % n, 2, u % y,.true.)    ! dU/dy
        call Grad_Mod_For_Phi(grid, u % n, 3, u % z,.true.)    ! dU/dz

        call Grad_Mod_For_Phi(grid, v % n, 1, v % x,.true.)    ! dV/dx
        call Grad_Mod_For_Phi(grid, v % n, 2, v % y,.true.)    ! dV/dy
        call Grad_Mod_For_Phi(grid, v % n, 3, v % z,.true.)    ! dV/dz

        call Grad_Mod_For_Phi(grid, w % n, 1, w % x,.true.)    ! dW/dx
        call Grad_Mod_For_Phi(grid, w % n, 2, w % y,.true.)    ! dW/dy
        call Grad_Mod_For_Phi(grid, w % n, 3, w % z,.true.)    ! dW/dz

        call Compute_Stresses(grid, dt, ini, uu)
        call Compute_Stresses(grid, dt, ini, vv)
        call Compute_Stresses(grid, dt, ini, ww)

        call Compute_Stresses(grid, dt, ini, uv)
        call Compute_Stresses(grid, dt, ini, uw)
        call Compute_Stresses(grid, dt, ini, vw)

        if(turbulence_model == REYNOLDS_STRESS_MODEL) then
          call Compute_F22(grid, ini, f22)
        end if

        call Compute_Stresses(grid, dt, ini, eps)

        call Calculate_Vis_T_Rsm(grid)

        if(heat_transfer == YES) then
          call Calculate_Heat_Flux(grid)
        end if
      end if

      if(turbulence_model == SPALART_ALLMARAS .or.  &
         turbulence_model == DES_SPALART) then
        call Calculate_Shear_And_Vorticity(grid)
        call Calculate_Vorticity(grid, u % n, v % n, w % n, vort)

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call Compute_Turbulent(grid, dt, ini, vis, n)
        call Calculate_Vis_T_Spalart_Allmaras(grid)
      end if

      ! Update the values at boundaries
      call Update_Boundary_Values(grid)

      ! End of the current iteration
      call Info_Mod_Iter_Print()

      if(ini >= min_ini) then
        if(coupling == 'SIMPLE') then
          call Control_Mod_Tolerance_For_Simple_Algorithm(simple_tol)
          if( u  % res <= simple_tol .and.  &
              v  % res <= simple_tol .and.  &
              w  % res <= simple_tol .and.  &
              mass_res <= simple_tol ) goto 4
        endif
      endif
    end do

    ! End of the current time step
4   call Info_Mod_Bulk_Print()

    ! Write the values in monitoring points
    do i=1,Nmon
      if(Cm(i)  > 0) then
        if(heat_transfer == NO) then
          write(10+i,'(I9,4E16.6)')                    &
           n, u % n(Cm(i)), v%n(Cm(i)), w%n(Cm(i)), p%n(Cm(i))
        else
          write(10+i,'(I9,5E16.6)')                    &
           n, u % n(Cm(i)), v%n(Cm(i)), w%n(Cm(i)), p%n(Cm(i)), t%n(Cm(i))
        end if
      end if
    end do

   if(PIPE==YES.or.JET==YES) then
     call CalcMn_Cylind(grid, n_stat, n)  !  calculate mean values
   else
     call Calculate_Mean(grid, n_stat, n)  !  calculate mean values
   end if

    !-----------------------------------------------------!
    !   Recalculate the pressure drop                     !
    !   to keep the constant mass flux                    !
    !                                                     !
    !   First Newtons law:                                !
    !                                                     !
    !   F = m * a                                         !
    !                                                     !
    !   where:                                            !
    !                                                     !
    !   a = dv / dt = dFlux / dt * 1 / (A * rho)          !
    !   m = rho * v                                       !
    !   F = Pdrop * l * A = Pdrop * v                     !
    !                                                     !
    !   finally:                                          !
    !                                                     !
    !   Pdrop * v = rho * v * dFlux / dt * 1 / (A * rho)  !
    !                                                     !
    !   after cancelling: v and rho, it yields:           !
    !                                                     !
    !   Pdrop = dFlux/dt/A                                !
    !-----------------------------------------------------!
    do m=1,grid % n_materials
      if( bulk(m) % flux_x_o /=  0.0 ) then
        bulk(m) % p_drop_x = (bulk(m) % flux_x_o - bulk(m) % flux_x)  &
                           / (dt * bulk(m) % area_x + TINY)
      end if
      if( bulk(m) % flux_y_o /=  0.0 ) then
        bulk(m) % p_drop_y = (bulk(m) % flux_y_o - bulk(m) % flux_y)  &
                           / (dt * bulk(m) % area_y + TINY)
      end if
      if( bulk(m) % flux_z_o /=  0.0 ) then
        bulk(m) % p_drop_z = (bulk(m) % flux_z_o - bulk(m) % flux_z)  &
                           / (dt * bulk(m) % area_z + TINY)
      end if
    end do

    !----------------------!
    !   Save the results   !
    !----------------------!
    inquire(file='exit_now', exist=exit_now)
    inquire(file='save_now', exist=save_now)

    ! Form the file name
    name_save = problem_name
    write(name_save(len_trim(problem_name)+1:                    &
                    len_trim(problem_name)+3), '(a3)'),   '-ts'
    write(name_save(len_trim(problem_name)+4:                    &
                    len_trim(problem_name)+9), '(i6.6)'), n

    ! Is it time to save the restart file?
    if(save_now .or. exit_now .or. mod(n,bsi) == 0) then
      call Comm_Mod_Wait
      call Save_Restart(grid, n, name_save)
    end if

    ! Is it time to save results for post-processing
    if(save_now .or. exit_now .or. mod(n,rsi) == 0) then
      call Comm_Mod_Wait
      call Save_Results(grid, name_save)

      ! Write results in user-customized format
      call User_Mod_Save_Results(grid, name_save)
    end if

    if(save_now) then
      open (9, file='save_now', status='old')
      close(9, status='delete')
    end if

    if(exit_now) then
      open (9, file='exit_now', status='old')
      close(9, status='delete')
    end if

  end do ! n, number of time steps

  if(this_proc < 2) then
    open(9, file='stop')
    close(9)
  end if

  if(this_proc  < 2) print *, '# Exiting !'

  !--------------------------------!
  !   Close the monitoring files   !
  !--------------------------------!
  do n = 1, Nmon
    if(Cm(n) > 0) close(10 + n)
  end do

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call Comm_Mod_End

  end program
