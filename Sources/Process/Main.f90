!==============================================================================!
  program Processor
!------------------------------------------------------------------------------!
!   Unstructured Finite Volume LES/RANS solver.                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Var_Mod
  use Solvers_Mod, only: D
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Correct_Velocity
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, m, n, Ndtt_temp
  real              :: Mres, CPUtim
  real              :: start, finish
  character(len=10) :: name_save
  logical           :: restar, multiple 
!---------------------------------[Interfaces]---------------------------------!
  interface
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine NewUVW(grid, var, Ui, dUidi, dUidj, dUidk,  &
                      Si, Sj, Sk, Di, Dj, Dk, Hi, dUjdi, dUkdi) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      integer         :: var
      type(Var_Type)  :: Ui
      real            :: dUidi(-grid % n_bnd_cells:grid % n_cells),  &
                         dUidj(-grid % n_bnd_cells:grid % n_cells),  &
                         dUidk(-grid % n_bnd_cells:grid % n_cells)
      real            :: Si(grid % n_faces),  &
                         Sj(grid % n_faces),  &
                         Sk(grid % n_faces) 
      real            :: Di(grid % n_faces),  &
                         Dj(grid % n_faces),  &
                         Dk(grid % n_faces) 
      real            :: Hi   (-grid % n_bnd_cells:grid % n_cells),  &
                         dUjdi(-grid % n_bnd_cells:grid % n_cells),  &
                         dUkdi(-grid % n_bnd_cells:grid % n_cells) 
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine Compute_Scalar(grid, var, phi, dphidx, dphidy, dphidz)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      integer         :: var
      type(Var_Type)  :: phi
      real            :: dphidx(-grid % n_bnd_cells:grid % n_cells),  &
                         dphidy(-grid % n_bnd_cells:grid % n_cells),  &
                         dphidz(-grid % n_bnd_cells:grid % n_cells)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine Compute_Turbulent(grid, var, phi, dphidx, dphidy, dphidz, Nstep)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      integer        :: var, Nstep
      type(Var_Type) :: phi
      real           :: dphidx(-grid % n_bnd_cells:grid % n_cells),  &
                        dphidy(-grid % n_bnd_cells:grid % n_cells),  &
                        dphidz(-grid % n_bnd_cells:grid % n_cells)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine Save_Gmv_Results(grid, namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      character, optional :: namAut*(*)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine Save_Dat_Results(grid, namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      character, optional :: namAut*(*)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
    subroutine Save_Restart(grid, namAut)  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 
      use all_mod
      use pro_mod
      use Grid_Mod
      implicit none
      type(Grid_Type)     :: grid
      character, optional :: namAut*(*)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine UserProbe2D(namAut)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      use all_mod
      use pro_mod
      implicit none
      character, optional :: namAut*(*)
    end subroutine
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine CalcShear(grid, Ui, Vi, Wi, She)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      use all_mod
      use pro_mod
      use les_mod
      use Grid_Mod
      implicit none
      type(Grid_Type) :: grid
      real            :: Ui(-grid % n_bnd_cells:grid % n_cells),  &
                         Vi(-grid % n_bnd_cells:grid % n_cells),  &
                         Wi(-grid % n_bnd_cells:grid % n_cells)
      real            :: She(-grid % n_bnd_cells:grid % n_cells)
    end subroutine
  end interface
!======================================================================!

  ! Grid used in computations
  type(Grid_Type) :: grid  

  call cpu_time(start)

  !--------------------------------!
  !   Splash out the logo screen   !             
  !--------------------------------!
  if(this_proc  < 2) then
    call logo
  endif

  !------------------------------!
  !   Start parallel execution   !
  !------------------------------!
  call StaPar()

  call Timex(CPUtim) 

  !--------------------------------------------!
  !   Open the command file, initialize line   !
  !    count to zero, and read problem name    ! 
  !--------------------------------------------!
  open(CMN_FILE, file='T-FlowS.cmn')    
  cmn_line_count = 0

  if(this_proc < 2) write(*,*) '# Input problem name:'
  call Tokenizer_Mod_Read_Line(CMN_FILE)  
  read(line % tokens(1), '(A80)')  name

  call Wait   

  ! Initialize parameters
  call IniPar
  
  ! Load the finite volume grid      
  call Load_Cns       (grid, this_proc)
  call Read_Problem           ! bad practice, should be avoided
  call Allocate_Memory(grid)
  call Load_Geo       (grid, this_proc)
  call BufLoa
  call Exchange(grid, grid % vol(-grid % n_bnd_cells))

  call wait   

  call Matrix_Mod_Topology(grid, A)
  call Matrix_Mod_Topology(grid, D)

  call Wait   

!>>>>>>>>>>>>>>>>>>>!
!                   !
!     TIME LOOP     !
!                   !
!<<<<<<<<<<<<<<<<<<<!
  Ndt  = 0
  Ndtt = 0
  call Load_Restart(grid, restar)

  if(restar) then
    call Load_Boundary_Conditions(grid, .FALSE.)
  end if

  ! Read command file (T-FlowS.cmn) 
  call ReaCom(grid, restar)

  ! Initialize variables
  if(.not. restar) then
    call Load_Boundary_Conditions(grid, .TRUE.)
    call Initialize_Variables(grid)
    call Wait
  end if   

  ! Interpolate between diff. meshes
  call Load_Ini(grid)

  ! Check if there are more materials
  multiple = .FALSE.
  i = StateMat(1)
  do m=1,grid % n_materials
    if(StateMat(m) /= i) multiple = .TRUE.
  end do

  if(this_proc  < 2)             &
    write(*,'(A18,15A11)')  &
    '#        N        ',   &
    ' Mass inb. ',          &
    '    U:     ',          &
    '    V:     ',          &
    '    W:     ',          &
    '    P:     '   

  if(this_proc  < 2)             &
    write(*,'(A18,15A11)')  &
    '  CFL max: ',          &
    '  Pe max:  ',          &
    '  Flux x:  ',          &
    '  Flux y:  ',          &
    '  Flux z:  ',          &
    'iterations ',          &
    '  dp/dx:   ',          &
    '    K:     '

  write(*,*) 'At line: ', cmn_line_count

  ! Loading data from previous computation   
  !  if(this_proc<2) write(*,*)'Reading data from previous computation on the same mesh'
  call Load_Restart_Ini(grid)

  ! Prepare ...
  call Compute_Geometry(grid)
  call Find_Bad        (grid)
  if(SIMULA==LES.and.MODE==SMAG.and..NOT.restar) call NearWallCell(grid)

  ! Prepare the gradient matrix for velocities
  call Compute_Gradient_Matrix(grid, .TRUE.) 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!        LET THE TIME        !
!     INTEGRATION  BEGIN     !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

  LinMon0 = ' '
  LinMon1 = '# [1]'
  LinMon2 = '# [2]'
  LinMon3 = '# [3]'
  LineRes = '#'

  if(ALGOR  ==  FRACT) call Pressure_Matrix_Fractional(grid)

  ! Print the areas of monitoring planes
  if(this_proc < 2) then
    do m=1,grid % n_materials
      write(*,'(A5,I2,A2,1PE12.3)') '# Ax(',m,')=', AreaX(m)
      write(*,'(A5,I2,A2,1PE12.3)') '# Ay(',m,')=', AreaY(m)
      write(*,'(A5,I2,A2,1PE12.3)') '# Az(',m,')=', AreaZ(m)
    end do
  end if

  if(Ndt == 0)  goto 6 

  do n=Ndtt+1, Ndt 

    Time = Time + dt

    if(Cm(1) > 0) then
      write(LinMon0(1: 6), '(I6)')      n;    
      write(LinMon0(7:18), '(1PE12.3)') Time; 
    end if

    if(SIMULA==DES_SPA) then
      call CalcShear(grid, U % n, V % n, W % n, Shear)
      call CalcVort (grid, U % n, V % n, W % n, Vort)
    end if

    if(SIMULA == LES) then
      call CalcShear(grid, U % n, V % n, W % n, Shear)
      if(MODE == DYN) call Compute_Sgs_Dynamic(grid) 
      if(MODE == WALE) call CalcWALE(grid) 
      call Compute_Sgs(grid)
    end if  

    If(SIMULA==HYB_ZETA) then  
      call Compute_Sgs_Dynamic(grid)      
      call Compute_Sgs_Hybrid(grid)
    end if

    call CalcConvect(grid)
    if(SIMULA==EBM.or.SIMULA==HJ) call CalcVISt_RSM(grid)
    
    do ini=1,Nini                   !  FRACTION & SIMPLE  

      if(.NOT. multiple) then 
        call GradP(grid, P % n, Px, Py, Pz)
      else 
        call GradP3(grid, P % n, Px, Py, Pz)
      end if
    
      call GraPhi(grid, U % n, 1, Ux, .TRUE.)    ! dU/dx
      call GraPhi(grid, U % n, 2, Uy, .TRUE.)    ! dU/dy
      call GraPhi(grid, U % n, 3, Uz, .TRUE.)    ! dU/dz
      call GraPhi(grid, V % n, 1, Vx, .TRUE.)    ! dV/dx
      call GraPhi(grid, V % n, 2, Vy, .TRUE.)    ! dV/dy
      call GraPhi(grid, V % n, 3, Vz, .TRUE.)    ! dV/dz
      call GraPhi(grid, W % n, 1, Wx, .TRUE.)    ! dW/dx
      call GraPhi(grid, W % n, 2, Wy, .TRUE.)    ! dW/dy 
      call GraPhi(grid, W % n, 3, Wz, .TRUE.)    ! dW/dz

      ! U velocity component
      call NewUVW(grid, 1,              &
                  U,                 &   
                  Ux,   Uy,   Uz,    & 
                  grid % sx,   grid % sy,   grid % sz,    &                             
                  grid % dx,   grid % dy,   grid % dz,    &                             
                  Px,   Vx,   Wx)      ! dP/dx, dV/dx, dW/dx

      ! V velocity component
      call NewUVW(grid, 2,              &
                  V,                 &
                  Vy,   Vx,   Vz,    & 
                  grid % sy,   grid % sx,   grid % sz,    &
                  grid % dy,   grid % dx,   grid % dz,    &
                  Py,   Uy,   Wy)      ! dP/dy, dU/dy, dW/dy

      ! W velocity component
      call NewUVW(grid, 3,              &
                  W,                 & 
                  Wz,   Wx,   Wy,    & 
                  grid % sz,   grid % sx,   grid % sy,    &                         
                  grid % dz,   grid % dx,   grid % dy,    &                         
                  Pz,   Uz,   Vz)      ! dP/dz, dU/dz, dV/dz

      if(ALGOR == FRACT) then
        call Exchange(grid, A % sav)  
        call Balance_Mass(grid)
        call Compute_Pressure_Fractional(grid)
      endif
      if(ALGOR == SIMPLE) then
        call Exchange(grid, A % sav)  
        call Balance_Mass(grid)
        call Compute_Pressure_Simple(grid)
      end if

      if(.NOT. multiple) then 
        call GradP(grid, PP % n,Px,Py,Pz)
      else 
        call GradP3(grid, PP % n,Px,Py,Pz)
      end if

      call Compute_Fluxes(grid)
      Mres = Correct_Velocity(grid) !  project the velocities

      ! Temperature
      if(HOT==YES) then
        call GraPhi(grid, T % n,1,phix,.TRUE.)            ! dT/dx
        call GraPhi(grid, T % n,2,phiy,.TRUE.)            ! dT/dy
        call GraPhi(grid, T % n,3,phiz,.TRUE.)            ! dT/dz
!!!!        call GraCorNew(grid, T % n,phix,phiy,phiz)    ! 2mat
        call Compute_Scalar(grid, 5, T,  &
                            phix, phiy, phiz)  ! dT/dx, dT/dy, dT/dz
      end if 

      ! Rans models
      if(SIMULA==K_EPS.or.SIMULA == HYB_PITM) then

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)
        call CalcShear(grid, U % n, V % n, W % n, Shear)
        call GraPhi(grid, Kin % n,1,phix,.TRUE.)             ! dK/dx
        call GraPhi(grid, Kin % n,2,phiy,.TRUE.)             ! dK/dy
        call GraPhi(grid, Kin % n,3,phiz,.TRUE.)             ! dK/dz
        call Compute_Turbulent(grid, 6, Kin,  &
                               phix, phiy, phiz, n) ! dK/dx, dK/dy, dK/dz

        call GraPhi(grid, Eps % n,1,phix,.TRUE.)           ! dEps/dx
        call GraPhi(grid, Eps % n,2,phiy,.TRUE.)           ! dEps/dy
        call GraPhi(grid, Eps % n,3,phiz,.TRUE.)           ! dEps/dz
        call Compute_Turbulent(grid, 7,  &
                               Eps, phix, phiy, phiz, n)
        call CalcVISt_KEps(grid)
      end if 

      if(SIMULA == K_EPS_VV .or.  &
         SIMULA == ZETA     .or.  &
         SIMULA == HYB_ZETA) then
        call CalcShear(grid, U % n, V % n, W % n, Shear)

        call GraPhi(grid, Kin % n,1,phix,.TRUE.)             ! dK/dx
        call GraPhi(grid, Kin % n,2,phiy,.TRUE.)             ! dK/dy
        call GraPhi(grid, Kin % n,3,phiz,.TRUE.)             ! dK/dz
        call Compute_Turbulent(grid, 6,   &
                               Kin, phix, phiy, phiz, n) ! dK/dx, dK/dy, dK/dz

        call GraPhi(grid, Eps % n,1,phix,.TRUE.)             ! dEps/dx
        call GraPhi(grid, Eps % n,2,phiy,.TRUE.)             ! dEps/dy
        call GraPhi(grid, Eps % n,3,phiz,.TRUE.)             ! dEps/dz
        call Compute_Turbulent(grid, 7, Eps, phix, phiy, phiz, n)
         
        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call GraPhi(grid, f22 % n,1,phix,.TRUE.)             ! df22/dx
        call GraPhi(grid, f22 % n,2,phiy,.TRUE.)             ! df22/dy
        call GraPhi(grid, f22 % n,3,phiz,.TRUE.)             ! df22/dz
        call Compute_F22(grid, 8, f22, phix, phiy, phiz, n) 

        call GraPhi(grid, v_2 % n,1,phix,.TRUE.)             ! dv_2/dx
        call GraPhi(grid, v_2 % n,2,phiy,.TRUE.)             ! dv_2/dy
        call GraPhi(grid, v_2 % n,3,phiz,.TRUE.)             ! dv_2/dz
        call Compute_Turbulent(grid, 9, v_2, phix, phiy, phiz, n)  

        call CalcVISt_KepsV2F(grid)
      end if                 

      if(SIMULA==EBM.or.SIMULA==HJ) then

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        if(SIMULA==EBM) call Time_And_Length_Scale(grid)

        call GraPhi(grid, U % n, 1, Ux,.TRUE.)    ! dU/dx
        call GraPhi(grid, U % n, 2, Uy,.TRUE.)    ! dU/dy
        call GraPhi(grid, U % n, 3, Uz,.TRUE.)    ! dU/dz
 
        call GraPhi(grid, V % n, 1, Vx,.TRUE.)    ! dV/dx
        call GraPhi(grid, V % n, 2, Vy,.TRUE.)    ! dV/dy
        call GraPhi(grid, V % n, 3, Vz,.TRUE.)    ! dV/dz

        call GraPhi(grid, W % n, 1, Wx,.TRUE.)    ! dW/dx
        call GraPhi(grid, W % n, 2, Wy,.TRUE.)    ! dW/dy 
        call GraPhi(grid, W % n, 3, Wz,.TRUE.)    ! dW/dz

        call GraPhi(grid, uu % n,1,phix,.TRUE.)             ! dK/dx
        call GraPhi(grid, uu % n,2,phiy,.TRUE.)             ! dK/dy
        call GraPhi(grid, uu % n,3,phiz,.TRUE.)             ! dK/dz
        call Compute_Stresses(grid, 6, uu, phix, phiy, phiz) ! dK/dx, dK/dy, dK/dz

        call GraPhi(grid, vv % n,1,phix,.TRUE.)             ! dEps/dx
        call GraPhi(grid, vv % n,2,phiy,.TRUE.)             ! dEps/dy
        call GraPhi(grid, vv % n,3,phiz,.TRUE.)             ! dEps/dz
        call Compute_Stresses(grid, 7, vv, phix, phiy, phiz)
         
        call GraPhi(grid, ww % n,1,phix,.TRUE.)             ! df22/dx
        call GraPhi(grid, ww % n,2,phiy,.TRUE.)             ! df22/dy
        call GraPhi(grid, ww % n,3,phiz,.TRUE.)             ! df22/dz
        call Compute_Stresses(grid, 8, ww, phix, phiy, phiz) 


        call GraPhi(grid, uv % n,1,phix,.TRUE.)             ! dv_2/dx
        call GraPhi(grid, uv % n,2,phiy,.TRUE.)             ! dv_2/dy
        call GraPhi(grid, uv % n,3,phiz,.TRUE.)             ! dv_2/dz
        call Compute_Stresses(grid, 9, uv, phix, phiy, phiz)  

        call GraPhi(grid, uw % n,1,phix,.TRUE.)             ! df22/dx
        call GraPhi(grid, uw % n,2,phiy,.TRUE.)             ! df22/dy
        call GraPhi(grid, uw % n,3,phiz,.TRUE.)             ! df22/dz
        call Compute_Stresses(grid, 10, uw, phix, phiy, phiz) 

        call GraPhi(grid, vw % n,1,phix,.TRUE.)             ! df22/dx
        call GraPhi(grid, vw % n,2,phiy,.TRUE.)             ! df22/dy
        call GraPhi(grid, vw % n,3,phiz,.TRUE.)             ! df22/dz
        call Compute_Stresses(grid, 11, vw, phix, phiy, phiz) 

        if(SIMULA==EBM) then
          call GraPhi(grid, f22 % n,1,phix,.TRUE.)             ! df22/dx
          call GraPhi(grid, f22 % n,2,phiy,.TRUE.)             ! df22/dy
          call GraPhi(grid, f22 % n,3,phiz,.TRUE.)             ! df22/dz
          call Compute_F22(grid, 12, f22, phix, phiy, phiz) 
        end if 

        call GraPhi(grid, Eps % n,1,phix,.TRUE.)             ! df22/dx
        call GraPhi(grid, Eps % n,2,phiy,.TRUE.)             ! df22/dy
        call GraPhi(grid, Eps % n,3,phiz,.TRUE.)             ! df22/dz
        call Compute_Stresses(grid, 13, Eps, phix, phiy, phiz) 
 
        call CalcVISt_RSM(grid)
      end if                 

      if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA) then
        call CalcShear(grid, U % n, V % n, W % n, Shear)
        call CalcVort(grid, U % n, V % n, W % n, Vort)

        ! Update the values at boundaries
        call Update_Boundary_Values(grid)

        call GraPhi(grid, VIS % n,1,phix,.TRUE.)             ! dVIS/dx
        call GraPhi(grid, VIS % n,2,phiy,.TRUE.)             ! dVIS/dy
        call GraPhi(grid, VIS % n,3,phiz,.TRUE.)             ! dVIS/dz
        call Compute_Turbulent(grid, 6,  &
                               VIS, phix, phiy, phiz, n)  ! dVIS/dx, dVIS/dy, dVIS/dz
        call CalcVISt_SPA_ALL(grid, n)
      end if

      ! Update the values at boundaries
      call Update_Boundary_Values(grid)

      write(LineRes(2:5),'(I4)') ini

      ! End of the current iteration 
      if(Cm(1) /= 0) write(*,'(A100)') LineRes     

      if(ALGOR == SIMPLE) then
        if( res(1) <= SIMTol .and. res(2) <= SIMTol .and. &
           res(3) <= SIMTol .and. res(4) <= SIMTol ) goto 4 
      endif
    end do 

    ! End of the current time step
4   if(Cm(1) /= 0) then
      write(*,'(A138)') LinMon0
      do m=1,grid % n_materials
        if(m .eq.1) write(*,'(A138)') LinMon1
        if(m .eq.2) write(*,'(A138)') LinMon2
        if(m .eq.3) write(*,'(A138)') LinMon3
      end do
    end if

    ! Write the values in monitoring points
    do i=1,Nmon
      if(Cm(i)  > 0) then                               
        if(HOT==NO) then
          write(10+i,'(I9,4E16.6)')                    &
           n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i))
        else
          write(10+i,'(I9,5E16.6)')                    &
           n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i)), T%n(Cm(i))
        end if
      end if
    end do  

   if(PIPE==YES.or.JET==YES) then 
     call CalcMn_Cylind(grid, Nstat, n)  !  calculate mean values 
!BOJAN     if(BUDG == YES.and.HOT==YES) call CalcBudgets_cylind(Nbudg, n)
!BOJAN     if(BUDG == YES.and.HOT==NO)  call CalcBudgets_cylind(Nbudg, n)
   else
     call Compute_Mean(grid, Nstat, n)  !  calculate mean values 
   end if

!-----------------------------------------------------!  
!                                                     ! 
!   Recalculate the pressure drop                     !
!   to keep the constant mass flux                    !
!                                                     !
!   First Newtons law:                                !
!   ~~~~~~~~~~~~~~~~~~                                !
!   F = m * a                                         !
!                                                     !
!   where:                                            !
!   ~~~~~~                                            !
!   a = dv / dt = dFlux / dt * 1 / (A * rho)          !
!   m = rho * V                                       !
!   F = Pdrop * l * A = Pdrop * V                     !
!                                                     !
!   finally:                                          !
!   ~~~~~~~~                                          !
!   Pdrop * V = rho * V * dFlux / dt * 1 / (A * rho)  !
!                                                     !
!   after cancelling: V and rho, it yields:           !
!                                                     !
!   Pdrop = dFlux/dt/A                                !
!                                                     !
!-----------------------------------------------------!
    do m=1,grid % n_materials
      if( FLUXoX(m)  /=  0.0 ) then
        PdropX(m) = (FLUXoX(m)-FLUXx(m)) / (dt*AreaX(m)+TINY) 
      end if
      if( FLUXoY(m)  /=  0.0 ) then
        PdropY(m) = (FLUXoY(m)-FLUXy(m)) / (dt*AreaY(m)+TINY) 
      end if
      if( FLUXoZ(m)  /=  0.0 ) then
        PdropZ(m) = (FLUXoZ(m)-FLUXz(m)) / (dt*AreaZ(m)+TINY) 
      end if
    end do

    ! Regular backup savings, each 1000 time steps            
    if(mod(n,1000) == 0) then                                
      Ndtt_temp = Ndtt
      Ndtt      = n
      name_save = 'savexxxxxx'
      write(name_save(5:10),'(I6.6)') n
      call Save_Restart(grid, name_save)                          
      Ndtt = Ndtt_temp
    end if   

    ! Is user forcing it to stop ?                    
    open(1, file='exit_now', err=7, status='old') 
      call Wait
      if(this_proc < 2) close(1, status='delete')
      Ndtt = n 
      goto 6 
7   continue 

    ! Is user forcing it to save ?                    
    open(1, file='save_now', err=8, status='old') 
      call Wait
      if(this_proc < 2) close(1, status='delete')
      Ndtt_temp = Ndtt
      Ndtt      = n
      name_save = 'savexxxxxx'
      write(name_save(5:10),'(I6.6)') n
      call Save_Restart    (grid, name_save)                          
      call Save_Gmv_Results(grid, name_save)
      call Save_Dat_Results(grid, name_save)
      call User_Save_Results(grid, n)  ! write results in user-customized format
      call SavParView(this_proc,grid % n_cells,name_save, Nstat, n)
      Ndtt = Ndtt_temp
8   continue 

  end do                    ! n, number of time steps  

  if(this_proc < 2) then
    open(9,file='stop')
    close(9)
  end if

  Ndtt = n - 1                                 

  !----------------------!
  !   Save the results   !
  !----------------------!
6 call Save_Restart     (grid)
  call Save_Ini         (grid)
  call Save_Gmv_Results (grid)     ! write results in GMV format. 
  call Save_Dat_Results (grid)     ! write results in FLUENT dat format. 
  call User_Save_Results(grid, n)  ! write results in user-customized format
  name_save = 'savexxxxxx'
  write(name_save(5:10),'(I6.6)') n
  call SavParView(this_proc, grid % n_cells, name_save, Nstat, n)
  call Wait

  if(this_proc  < 2) write(*,*) '# Exiting !'

  !----------------------------!
  !   Close the command file   !
  !----------------------------!
  close(CMN_FILE)                

  !--------------------------------!
  !   Close the monitoring files   !
  !--------------------------------!
  do n=1,Nmon
    if(Cm(n)  > 0) close(10+n)
  end do

  !----------------------------!
  !   End parallel execution   !
  !----------------------------!
  call endpar            
  call cpu_time(finish)
  if(this_proc < 2) print '("Time = ",f14.3," seconds.")',finish-start

  end program
