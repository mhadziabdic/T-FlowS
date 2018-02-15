!==============================================================================!
  subroutine Compute_F22(grid, ini, var, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves eliptic relaxation equations for f22.               !
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
  integer         :: ini
  integer         :: var
  type(Var_Type)  :: phi
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter, miter
  real              :: Fex, Fim 
  real              :: A0, A12, A21
  real              :: error, tol
  real              :: phi_x_f, phi_y_f, phi_z_f
  character(len=80) :: coupling
  character(len=80) :: precond
  character(len=80) :: td_inertia    ! time-disretization for inerita  
  character(len=80) :: td_advection  ! time-disretization for advection
  character(len=80) :: td_diffusion  ! time-disretization for diffusion 
  character(len=80) :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor                 
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
!==============================================================================! 
!
!   The form of equations which are solved:
!
!      /           /              /
!     | df22      | f22 dV       | f22hg dV
!   - | ---- dS + | ------   =   | --------
!     |  dy       |  Lsc^2       |  Lsc^2
!    /           /              /
!
!   Dimension of the system under consideration
!
!     [A]{f22} = {b}   [kg K/s]
!
!   Dimensions of certain variables:
!
!     f22            [1/s]
!     Lsc            [m]
!
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

  ! Old values (o) and older than old (oo)
  if(ini == 1) then
    do c=1,grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
      phi % d_oo(c) = phi % d_o(c)
      phi % d_o (c) = 0.0 
      phi % c_oo(c) = phi % c_o(c)
      phi % c_o (c) = phi % c(c) 
    end do
  end if

  ! New values
  do c=1,grid % n_cells
    phi % c(c) = 0.0
  end do

  ! Gradients
  call GraPhi(grid, phi % n, 1, phi_x, .true.)
  call GraPhi(grid, phi % n, 2, phi_y, .true.)
  call GraPhi(grid, phi % n, 3, phi_z, .true.)

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

    phi_x_f = fF(s) * phi_x(c1) + (1.0 - fF(s)) * phi_x(c2)
    phi_y_f = fF(s) * phi_y(c1) + (1.0 - fF(s)) * phi_y(c2)
    phi_z_f = fF(s) * phi_z(c1) + (1.0 - fF(s)) * phi_z(c2)


    ! Total (exact) diffusive flux
    Fex = (  phi_x_f * grid % sx(s)   &
           + phi_y_f * grid % sy(s)   &
           + phi_z_f * grid % sz(s) )

    A0 = Scoef(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: Scoef is
    !  not corrected at interface between materials)
    Fim=(   phi_x_f * grid % dx(s)        &
          + phi_y_f * grid % dy(s)        &
          + phi_z_f * grid % dz(s)) * A0

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

    ! Calculate the coefficients for the sysytem matrix
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

      ! Fill the system matrix
      if(c2  > 0) then
        A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A12
        A % val(A % dia(c1))  = A % val(A % dia(c1))  + A12
        A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A21
        A % val(A % dia(c2))  = A % val(A % dia(c2))  + A21
      else if(c2  < 0) then
        ! Outflow is not included because it was causing problems  
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW)) then                    
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1) = b(c1) + A12 * phi % n(c2)

        else

        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL).or.                          &
            (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) ) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          !---------------------------------------------------------------!
          !   Source coefficient is filled in SourceF22.f90 in order to   !
          !   get updated values of f22 on the wall.  Otherwise f22       !
          !   equation does not converge very well                        !
          !   b(c1) = b(c1) + A12 * phi % n(c2)                           !
          !---------------------------------------------------------------!
        else if( Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) then  
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          A % bou(c2) = - A12  ! cool parallel stuff
        endif
      end if     
     end if
    end if

  end do  ! through faces

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion == 'ADAMS_BASHFORTH') then 
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * phi % d_o(c) - 0.5 * phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion == 'CRANK_NICOLSON') then 
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * phi % d_o(c)
    end do  
  end if
                 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff == 'ADAMS_BASHFORTH') then
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * phi % c_o(c) - 0.5 * phi % c_oo(c)
    end do 
  end if

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff == 'CRANK_NICOLSON') then
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * phi % c(c) + 0.5 * phi % c_o(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff == 'FULLY_IMPLICIT') then
    do c=1,grid % n_cells
      b(c) = b(c) + phi % c(c)
    end do 
  end if

  !-------------------------------------!
  !                                     !  
  !   Source terms and wall function    !
  !   (Check if it is good to call it   !
  !    before the under relaxation ?)   !
  !                                     !
  !-------------------------------------!
  if(turbulence_model == 'EBM') then
    call Source_F22_Ebm(grid)
  else
    call Source_F22_K_Eps_V2_F(grid)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!

  ! Type of coupling is important
  call Control_Mod_Pressure_Momentum_Coupling(coupling)

  ! Set under-relaxation factor
  urf = 1.0
  if(coupling == 'SIMPLE')  &
    call Control_Mod_Simple_Underrelaxation_For_Turbulence(urf)

  do c=1,grid % n_cells
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - urf)*phi % n(c) / urf
    A % val(A % dia(c)) = A % val(A % dia(c)) / urf
  end do 

  ! Get tolerance for linear solver
  call Control_Mod_Tolerance_For_Turbulence_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of solver iterations on coupling method
  if(coupling == 'PROJECTION') miter = 300
  if(coupling == 'SIMPLE')     miter =   5

  niter=miter
  call cg(A, phi % n, b, precond, niter, tol, res(var), error)
  
  call Info_Mod_Iter_Fill_At(3, 4, phi % name, niter, res(var))

  call Exchange(grid, phi % n)

  end subroutine
