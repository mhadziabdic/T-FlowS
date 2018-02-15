!==============================================================================!
  subroutine Compute_Turbulent(grid, dt, ini, var, phi, Nstep)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equations for different turbulent         !
!   variables.                                                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
  use Var_Mod
  use Grid_Mod
  use Info_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use Work_Mod,    only: phi_x => r_cell_01,  &
                         phi_y => r_cell_02,  &
                         phi_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
  integer         :: ini
  integer         :: var
  type(Var_Type)  :: phi
  integer         :: Nstep
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter, miter, mat
  real              :: Fex, Fim 
  real              :: phis
  real              :: A0, A12, A21
  real              :: error, tol
  real              :: vis_eff
  real              :: phi_x_f, phi_y_f, phi_z_f
  character(len=80) :: coupling
  character(len=80) :: precond
  character(len=80) :: sd_advection  ! advection scheme
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  character(len=80) :: td_inertia    ! time-disretization for inerita  
  character(len=80) :: td_advection  ! time-disretization for advection
  character(len=80) :: td_diffusion  ! time-disretization for diffusion 
  character(len=80) :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor                 
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
!==============================================================================!
!                                                                              ! 
!   The form of equations which are solved:                                    !   
!                                                                              !
!      /               /                /                     /                !
!     |     dphi      |                | mu_eff              |                 !
!     | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV            !
!     |      dt       |                |  sigma              |                 !
!    /               /                /                     /                  !
!                                                                              !
!------------------------------------------------------------------------------!

  A % val = 0.0

  b=0.0

  ! This is important for "copy" boundary conditions. Find out why !
  do c=-grid % n_bnd_cells,-1
    A % bou(c)=0.0
  end do

  !-------------------------------------! 
  !   Initialize variables and fluxes   !
  !-------------------------------------! 

  call Control_Mod_Time_Integration_For_Inertia(td_inertia)
  call Control_Mod_Time_Integration_For_Advection(td_advection)
  call Control_Mod_Time_Integration_For_Diffusion(td_diffusion)
  call Control_Mod_Time_Integration_For_Cross_Diffusion(td_cross_diff)

  call Control_Mod_Turbulence_Model(turbulence_model)
  call Control_Mod_Turbulence_Model_Variant(turbulence_model_variant)

  ! Old values (o) and older than old (oo)
  if(ini == 1) then
    do c = 1, grid % n_cells
      phi % oo(c)  = phi % o(c)
      phi % o (c)  = phi % n(c)
      phi % a_oo(c) = phi % a_o(c)
      phi % a_o (c) = 0.0 
      phi % d_oo(c) = phi % d_o(c)
      phi % d_o (c) = 0.0 
      phi % c_oo(c) = phi % c_o(c)
      phi % c_o (c) = phi % c(c) 
    end do
  end if

  ! Gradients
  call GraPhi(grid, phi % n, 1, phi_x, .true.)
  call GraPhi(grid, phi % n, 2, phi_y, .true.)
  call GraPhi(grid, phi % n, 3, phi_z, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Turbulence(sd_advection)
  call Control_Mod_Blending_Coefficient_Turbulence(blend)
  
  ! Compute phimax and phimin
  do mat = 1, grid % n_materials
    if(sd_advection .ne. 'CENTRAL') then
      call Compute_Minimum_Maximum(grid, phi % n)  ! or phi % o ???
      goto 1
    end if
  end do

  ! New values
1 do c = 1, grid % n_cells
    phi % a(c)    = 0.0
    phi % c(c)    = 0.0
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s) 

    ! Velocities on "orthogonal" cell centers 
    if(c2 > 0 .or.  &
       c2 < 0.and.Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then
      phis =        grid % f(s)  * phi % n(c1)   &
           + (1.0 - grid % f(s)) * phi % n(c2)

      ! Compute phis with desired advection scheme
      if(sd_advection .ne. 'CENTRAL') then
        call Advection_Scheme(grid, phis, s, phi % n,           &
                              phi_x, phi_y, phi_z,              &
                              grid % dx, grid % dy, grid % dz,  &
                              sd_advection, blend) 
      end if 

      ! Compute advection term
      if(ini == 1) then 
        if(c2  > 0) then
          phi % a_o(c1)=phi % a_o(c1) - Flux(s) * phis
          phi % a_o(c2)=phi % a_o(c2) + Flux(s) * phis
        else
          phi % a_o(c1)=phi % a_o(c1) - Flux(s) * phis
        endif 
      end if
      if(c2  > 0) then
        phi % a(c1)=phi % a(c1)-Flux(s) * phis
        phi % a(c2)=phi % a(c2)+Flux(s) * phis
      else
        phi % a(c1)=phi % a(c1)-Flux(s) * phis
      endif 

      ! Store upwinded part of the advection term in "c"
      if(coupling .ne. 'PROJECTION') then
        if(Flux(s)  < 0) then   ! from c2 to c1
        phi % c(c1) = phi % c(c1) - Flux(s) * phi % n(c2)
          if(c2  > 0) then
            phi % c(c2) = phi % c(c2) + Flux(s) * phi % n(c2)
          endif
        else 
          phi % c(c1) = phi % c(c1) - Flux(s) * phi % n(c1)
          if(c2  > 0) then
            phi % c(c2) = phi % c(c2) + Flux(s) * phi % n(c1)
          endif
        end if
      end if   
    end if     ! c2 > 0 
  end do    ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(td_advection == 'ADAMS_BASHFORTH') then
    do c = 1, grid % n_cells
      b(c) = b(c) + (1.5*phi % a_o(c) - 0.5*phi % a_oo(c) - phi % c(c))
    end do  
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(td_advection == 'CRANK_NICOLSON') then
    do c = 1, grid % n_cells
      b(c) = b(c) + (0.5 * ( phi % a(c) + phi % a_o(c) ) - phi % c(c))
    end do  
  endif

  ! Fully implicit treatment of convective fluxes 
  if(td_advection == 'FULLY_IMPLICIT') then
    do c = 1, grid % n_cells
      b(c) = b(c) + (phi % a(c) - phi % c(c))
    end do  
  end if     
          
  ! New values
  do c = 1, grid % n_cells
    phi % c(c) = 0.0
  end do
  
  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces       

    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s)   

    vis_eff = VISc + (fF(s)*vis_t(c1) + (1.0-fF(s))*vis_t(c2))/phi % Sigma 

    if(turbulence_model == 'SPA_ALL' .or.  &
       turbulence_model == 'DES_SPA')      &
      vis_eff = VISc+(fF(s)*VIS % n(c1)+(1.0-fF(s))*VIS % n(c2))  &
             / phi % Sigma

    if(turbulence_model == 'HYB_ZETA')                                    &
      vis_eff = VISc + (fF(s)*vis_t_eff(c1) + (1.0-fF(s))*vis_t_eff(c2))  &
             / phi % Sigma

    phi_x_f = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2)
    phi_y_f = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
    phi_z_f = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)

    ! This implements zero gradient for k
    if(turbulence_model == 'K_EPS' .and.  &
       turbulence_model_variant == 'HIGH_RE') then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then  
          phi_x_f = 0.0 
          phi_y_f = 0.0
          phi_z_f = 0.0
          vis_eff  = 0.0
        end if 
      end if
    end if

    if(turbulence_model == 'ZETA'     .or.  &
       turbulence_model == 'K_EPS_VV' .or.  &
       turbulence_model == 'HYB_ZETA') then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          if(sqrt(TauWall(c1))*grid % wall_dist(c1)/VISc>2.0) then      
            phi_x_f = 0.0
            phi_y_f = 0.0
            phi_z_f = 0.0
            vis_eff  = 0.0
          end if
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    Fex = vis_eff * (  phi_x_f * grid % sx(s)  &
                    + phi_y_f * grid % sy(s)  &
                    + phi_z_f * grid % sz(s) )

    A0 = vis_eff * Scoef(s) 

    ! Implicit diffusive flux
    ! (this is a very crude approximation: Scoef is
    !  not corrected at interface between materials)
    Fim = (  phi_x_f * grid % dx(s)                      &
           + phi_y_f * grid % dy(s)                      &
           + phi_z_f * grid % dz(s) ) * A0

    ! Straight diffusion part 
    if(ini == 1) then
      if(c2  > 0) then
        phi % d_o(c1) = phi % d_o(c1) + (phi % n(c2)-phi % n(c1))*A0   
        phi % d_o(c2) = phi % d_o(c2) - (phi % n(c2)-phi % n(c1))*A0    
      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) /= SYMMETRY) then
          phi % d_o(c1) = phi % d_o(c1) + (phi % n(c2)-phi % n(c1))*A0   
        end if 
      end if 
    end if

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + Fex - Fim 
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - Fex + Fim 
    end if 

    ! Compute coefficients for the sysytem matrix
    if( (td_diffusion == 'CRANK_NICOLSON') .or.  &
        (td_diffusion == 'FULLY_IMPLICIT') ) then  

      if(td_diffusion == 'CRANK_NICOLSON') then  
        A12 = 0.5 * A0 
        A21 = 0.5 * A0 
      end if

      if(td_diffusion == 'FULLY_IMPLICIT') then 
        A12 = A0 
        A21 = A0
      end if

      if(coupling .ne. 'PROJECTION') then
        A12 = A12  - min(Flux(s), real(0.0)) 
        A21 = A21  + max(Flux(s), real(0.0))
      endif

      ! Fill the system matrix
      if(c2  > 0) then
        A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A12
        A % val(A % dia(c1))  = A % val(A % dia(c1))  + A12
        A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A21
        A % val(A % dia(c2))  = A % val(A % dia(c2))  + A21
      else if(c2  < 0) then

        ! Outflow is not included because it was causing problems     
        if((Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW)  .or.                 &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL)    .or.                 &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE).or.                 &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT) .or.                 &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) ) then                               
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1) = b(c1) + A12 * phi % n(c2)
        else if( Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) then  
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          A % bou(c2) = - A12  ! cool parallel stuff
        endif
      end if     
    end if
  end do  ! through faces

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion == 'ADAMS_BASHFORTH') then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % d_o(c) - 0.5 * phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion == 'CRANK_NICOLSON') then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % d_o(c)
    end do  
  end if
             
  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff == 'ADAMS_BASHFORTH') then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % c_o(c) - 0.5 * phi % c_oo(c)
    end do 
  end if
    
  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff == 'CRANK_NICOLSON') then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % c(c) + 0.5 * phi % c_o(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff == 'FULLY_IMPLICIT') then
    do c = 1, grid % n_cells
      b(c) = b(c) + phi % c(c)
    end do 
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; linear interpolation
  if(td_inertia == 'LINEAR') then
    do c = 1, grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia == 'PARABOLIC') then
    do c = 1, grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * A0
      b(c) = b(c) + 2.0*A0 * phi % o(c) - 0.5*A0 * phi % oo(c)
    end do
  end if

  !-------------------------------------!
  !                                     !  
  !   Source terms and wall function    !
  !   (Check if it is good to call it   !
  !    before the under relaxation ?)   !
  !                                     !
  !-------------------------------------!
  if(turbulence_model == 'K_EPS') then 
    if(phi % name == 'KIN') call Source_Kin_K_Eps(grid)
    if(phi % name == 'EPS') call Source_Eps_K_Eps(grid)
  end if

  if(turbulence_model == 'K_EPS_VV' .or.  &
     turbulence_model == 'ZETA'     .or.  &
     turbulence_model == 'HYB_ZETA') then
    if(phi % name == 'KIN') call Source_Kin_K_Eps_V2_F(grid)
    if(phi % name == 'EPS') call Source_Eps_K_Eps_V2_F(grid)
    if(phi % name == 'V^2') call Source_V2_K_Eps_V2_F(grid, Nstep)
  end if

  if(turbulence_model == 'SPA_ALL' .or.  &
     turbulence_model == 'DES_SPA')      &
    call Source_Vis_Spalart_Almaras(grid, phi_x, phi_y, phi_z)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!

  ! Type of coupling is important
  call Control_Mod_Pressure_Momentum_Coupling(coupling)

  ! Set under-relaxation factor
  urf = 1.0
  if(coupling == 'SIMPLE')   &
    call Control_Mod_Simple_Underrelaxation_For_Turbulence(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - urf) * phi % n(c) / urf
    A % val(A % dia(c)) = A % val(A % dia(c)) / urf
  end do

  ! Get tolerance for linear solvers
  call Control_Mod_Tolerance_For_Turbulence_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  if(coupling == 'PROJECTION') miter = 10 
  if(coupling == 'SIMPLE')     miter =  5

  niter=miter
  call cg(A, phi % n, b, precond, niter, tol, res(var), error)

  do c = 1, grid % n_cells
    if( phi%n(c)<0.0)then
      phi % n(c) = phi % o(c)
    end if
  end do 

  if(turbulence_model == 'K_EPS'    .or.  &
     turbulence_model == 'K_EPS_VV' .or.  &
     turbulence_model == 'ZETA'     .or.  &
     turbulence_model == 'HYB_ZETA') then
    if(phi % name == 'KIN')  &
      call Info_Mod_Iter_Fill_At(3, 1, phi % name, niter, res(var))
    if(phi % name == 'EPS')  &
      call Info_Mod_Iter_Fill_At(3, 2, phi % name, niter, res(var))
    if(phi % name == 'V^2')  &
      call Info_Mod_Iter_Fill_At(3, 3, phi % name, niter, res(var))
  end if

  call Exchange(grid, phi % n)

  end subroutine
