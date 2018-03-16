!==============================================================================!
  subroutine Compute_Temperature(grid, dt, ini, phi)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use rans_mod
  use Comm_Mod
  use Var_Mod
  use Grid_Mod
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use Work_Mod,    only: phi_x       => r_cell_01,  &
                         phi_y       => r_cell_02,  &
                         phi_z       => r_cell_03,  &
                         u1uj_phij   => r_cell_04,  &
                         u2uj_phij   => r_cell_05,  &
                         u3uj_phij   => r_cell_06,  &
                         u1uj_phij_x => r_cell_07,  &
                         u2uj_phij_y => r_cell_08,  &
                         u3uj_phij_z => r_cell_09    
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Grid_Type) :: grid
  integer         :: ini
  real            :: dt
  type(Var_Type)  :: phi
!----------------------------------[Calling]-----------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------! 
  integer           :: n, c, s, c1, c2, niter, mat
  real              :: A0, A12, A21
  real              :: ini_res, tol
  real              :: CONeff1, FUex1, FUim1, phixS1, phiyS1, phizS1
  real              :: CONeff2, FUex2, FUim2, phixS2, phiyS2, phizS2
  real              :: Stot, phis, Prt, Prt1, Prt2
  character(len=80) :: coupling   ! pressure-momentum coupling
  character(len=80) :: precond    ! preconditioner
  integer           :: adv_scheme  ! space-discretiztion of advection scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  integer           :: td_inertia    ! time-disretization for inerita  
  integer           :: td_advection  ! time-disretization for advection
  integer           :: td_diffusion  ! time-disretization for diffusion 
  integer           :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor                 
!------------------------------------------------------------------------------!
!     
!  The form of equations which are solved:    
!     
!     /                /                 /            
!    |        dT      |                 |             
!    | rho Cp -- dV   | rho u Cp T dS = | lambda  DIV T dS 
!    |        dt      |                 |             
!   /                /                 /              
!
!
!  Dimension of the system under consideration
!
!     [A]{T} = {b}   [J/s = W]  
!
!  Dimensions of certain variables:
!
!     Cp     [J/kg K]
!     lambda [W/m K]
!
!     A      [kg/s]
!     T      [K]
!     b      [kg K/s] 
!     Flux   [kg/s]
!     CT*,   [kg K/s] 
!     DT*,   [kg K/s] 
!     XT*,   [kg K/s]
! 
!==============================================================================!

!  if(turbulence_model == REYNOLDS_STRESS_MODEL) then
!    TDC = 1.0         
!  else
!    TDC = 1.0       
!  end if        

  do n = 1, A % row(grid % n_cells+1) ! to je broj nonzero + 1
    A % val(n) = 0.0
  end do
  A % val = 0.0

  do c = 1, grid % n_cells
    b(c)=0.0
  end do

  ! This is important for "copy" boundary conditions. Find out why !
  ! (I just coppied this from NewUVW.f90. I don't have a clue if it
  ! is of any importance at all. Anyway, I presume it can do no harm.)
  do c = -grid % n_bnd_cells,-1
    A % bou(c)=0.0
  end do

  !-------------------------------------! 
  !   Initialize variables and fluxes   !
  !-------------------------------------! 

  call Control_Mod_Time_Integration_For_Inertia(td_inertia)
  call Control_Mod_Time_Integration_For_Advection(td_advection)
  call Control_Mod_Time_Integration_For_Diffusion(td_diffusion)
  call Control_Mod_Time_Integration_For_Cross_Diffusion(td_cross_diff)

  ! Old values (o and oo)
  if(ini.lt.2) then
    do c = 1, grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
      phi % a_oo(c) = phi % a_o(c)
      phi % a_o (c) = 0.0
      phi % d_oo(c) = phi % d_o(c)
      phi % d_o (c) = 0.0 
      phi % c_oo(c) = phi % c_o(c)
      phi % c_o (c) = phi % c(c) 
    end do
  end if

  ! Gradients
  call Grad_Mod_For_Phi(grid, phi % n, 1, phi_x, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 2, phi_y, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 3, phi_z, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Energy(adv_scheme)
  call Control_Mod_Blending_Coefficient_Energy(blend)
  
  ! Compute phimax and phimin
  do mat = 1, grid % n_materials
    if(adv_scheme .ne. CENTRAL) then
      call Calculate_Minimum_Maximum(grid, phi % n)  ! or phi % o ???
      goto 1  ! why this???
    end if
  end do

  ! New values
1 do c = 1, grid % n_cells
    phi % a(c) = 0.0
    phi % c(c) = 0.0  ! use phi % c for upwind advective fluxes
  end do

  !----------------------------------!
  !   Browse through all the faces   !
  !----------------------------------!
  do s=1,grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    phis =      grid % f(s)  * phi % n(c1)   &
         + (1.0-grid % f(s)) * phi % n(c2)

    ! Compute phis with desired advection scheme
    if(adv_scheme .ne. CENTRAL) then
      call Advection_Scheme(grid, phis, s, phi % n,           &
                            phi_x, phi_y, phi_z,              &
                            grid % dx, grid % dy, grid % dz,  &
                            adv_scheme, blend) 
    end if

    ! Compute advection term
    if(ini.eq.1) then
      if(c2.gt.0) then
        phi % a_o(c1) = phi % a_o(c1)-flux(s)*phis*capacity
        phi % a_o(c2) = phi % a_o(c2)+flux(s)*phis*capacity
      else
        phi % a_o(c1) = phi % a_o(c1)-flux(s)*phis*capacity
      endif
    end if
    if(c2.gt.0) then
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
      phi % a(c2) = phi % a(c2)+flux(s)*phis*capacity
    else
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
    endif

    ! Store upwinded part of the advection term in "c"
    if(coupling .ne. 'PROJECTION') then
      if(flux(s).lt.0) then   ! from c2 to c1
        phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c2) * capacity
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c2) * capacity
        endif
      else
        phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c1) * capacity
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c1) * capacity
        endif
      end if
    end if  

 
  end do  ! through sides

  !-----------------------------!
  !                             !
  !   Temporal discretization   !
  !                             !
  !-----------------------------!
   
  ! Adams-Bashforth scheeme for advection fluxes
  if(td_advection == ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5*phi % a_o(c) - 0.5*phi % a_oo(c) - phi % c(c)
    end do
  endif

  ! Crank-Nicholson scheeme for advection fluxes
  if(td_advection == CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * (phi % a(c) + phi % a_o(c)) - phi % c(c)
    end do
  endif

  ! Fully implicit treatment of advection fluxes
  if(td_advection == FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + phi % a(c) - phi % c(c)
    end do
  end if

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  ! Set phi % c back to zero 
  do c = 1, grid % n_cells
    phi % c(c) = 0.0  
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  Prt = 0.9

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
     
    if(turbulence_model /= LES .or.  &
       turbulence_model /= DNS) then
      Prt1 = 1.0/( 0.5882 + 0.228*(vis_t(c1)/(viscosity+1.0e-12)) - 0.0441*                  &
            (vis_t(c1)/(viscosity+1.0e-12))**2.0*(1.0 - exp(-5.165*( viscosity/(vis_t(c1)+1.0e-12) ))) )
      Prt2 = 1.0/( 0.5882 + 0.228*(vis_t(c2)/(viscosity+1.0e-12)) - 0.0441*                  &
           (vis_t(c2)/(viscosity+1.0e-12))**2.0*(1.0 - exp(-5.165*( viscosity/(vis_t(c2)+1.0e-12) ))) )
      Prt = fw(s)*Prt1 + (1.0-fw(s))*Prt2
    end if

    ! Gradients on the cell face 
    if(c2 > 0 .or.  &
       c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then
      if(grid % material(c1) == grid % material(c2)) then
        phixS1 = fw(s)*phi_x(c1) + (1.0-fw(s))*phi_x(c2) 
        phiyS1 = fw(s)*phi_y(c1) + (1.0-fw(s))*phi_y(c2)
        phizS1 = fw(s)*phi_z(c1) + (1.0-fw(s))*phi_z(c2)
        phixS2 = phixS1 
        phiyS2 = phiyS1 
        phizS2 = phizS1 
        CONeff1 =      grid % f(s) * ( conductivity                 &
                                     + capacity*vis_t(c1)/Prt )  &
                + (1.-grid % f(s)) * ( conductivity                 &
                                     + capacity*vis_t(c2)/Prt )
        CONeff2 = CONeff1 
      else 
        phixS1 = phi_x(c1) 
        phiyS1 = phi_y(c1) 
        phizS1 = phi_z(c1) 
        phixS2 = phi_x(c2) 
        phiyS2 = phi_y(c2) 
        phizS2 = phi_z(c2) 
        CONeff1 =   conductivity                 &
                  + capacity*vis_t(c1)/Prt   
        CONeff2 =   conductivity                 &
                  + capacity*vis_t(c2)/Prt   
      end if
    else
      phixS1 = phi_x(c1) 
      phiyS1 = phi_y(c1) 
      phizS1 = phi_z(c1) 
      phixS2 = phixS1 
      phiyS2 = phiyS1 
      phizS2 = phizS1 
      CONeff1 =   conductivity                 &
                + capacity*vis_t(c1)/Prt   
      CONeff2 = CONeff1 
    endif


    if(turbulence_model == K_EPS_ZETA_F  .or.  &
       turbulence_model == K_EPS         .or.  &
       turbulence_model == HYBRID_K_EPS_ZETA_F) then
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          CONeff1 = con_wall(c1)
          CONeff2 = CONeff1
        end if
      end if
    end if  

    ! Total (exact) diffusive flux
    FUex1 = CONeff1*(  phixS1*grid % sx(s)  &
                     + phiyS1*grid % sy(s)  &
                     + phizS1*grid % sz(s))
    FUex2 = CONeff2*(  phixS2*grid % sx(s)  &
                     + phiyS2*grid % sy(s)  &
                     + phizS2*grid % sz(s))

    ! Implicit diffusive flux
    FUim1 = CONeff1*f_coef(s)*    &
            (  phixS1*grid % dx(s)      &
             + phiyS1*grid % dy(s)      &
             + phizS1*grid % dz(s) )
    FUim2 = CONeff2*f_coef(s)*    &
            (  phixS2*grid % dx(s)      &
             + phiyS2*grid % dy(s)      &
             + phizS2*grid % dz(s) )

    ! Straight diffusion part 
    if(ini.lt.2) then
      if(c2.gt.0) then
        if(grid % material(c1) == grid % material(c2)) then
          phi % d_o(c1) = phi % d_o(c1)  &
                       + CONeff1*f_coef(s)*(phi % n(c2) - phi % n(c1)) 
          phi % d_o(c2) = phi % d_o(c2)  &
                       - CONeff2*f_coef(s)*(phi % n(c2) - phi % n(c1))   
        else
          phi % d_o(c1) = phi % d_o(c1)            &
                       + 2.*conductivity   &
                       * f_coef(s) * (phi_face(s) - phi % n(c1)) 
          phi % d_o(c2) = phi % d_o(c2)            &
                       - 2.*conductivity   &
                       * f_coef(s)*(phi % n(c2) - phi_face(s))   
        end if
      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2).ne.SYMMETRY) then 
          if(grid % material(c1) == grid % material(c2)) then
            phi % d_o(c1) = phi % d_o(c1)  &
                + CONeff1*f_coef(s)*(phi % n(c2) - phi % n(c1))   
          else
            phi % d_o(c1) = phi % d_o(c1)  &
                + 2.*conductivity*f_coef(s)*(phi_face(s)-phi % n(c1)) 
          end if
        end if
      end if 
    end if

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + FUex1 - FUim1 
    if(c2.gt.0) then
      phi % c(c2) = phi % c(c2) - FUex2 + FUim2 
    end if 

    ! Calculate the coefficients for the sysytem matrix
    if( (td_diffusion == CRANK_NICOLSON) .or.  &
        (td_diffusion == FULLY_IMPLICIT) ) then

      if(td_diffusion == CRANK_NICOLSON) then 
        if(grid % material(c1) == grid % material(c2)) then
          A12 = .5*CONeff1*f_coef(s)  
          A21 = .5*CONeff2*f_coef(s)  
        else
          A12 = conductivity*f_coef(s)  
          A21 = conductivity*f_coef(s)  
        end if
      end if

      if(td_diffusion == FULLY_IMPLICIT) then 
        if(grid % material(c1) == grid % material(c2)) then
          A12 = CONeff1*f_coef(s)  
          A21 = CONeff2*f_coef(s)  
        else
          A12 = 2.*conductivity*f_coef(s)  
          A21 = 2.*conductivity*f_coef(s)  
        end if
      end if

      if(coupling .ne. 'PROJECTION') then
        A12 = A12  - min(flux(s), 0.0) * capacity
        A21 = A21  + max(flux(s), 0.0) * capacity
      endif
                
      ! Fill the system matrix
      if(c2.gt.0) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + A21
        if(grid % material(c1) == grid % material(c2)) then
          A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A12
          A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A21
        else
          b(c1) = b(c1) + A12*phi_face(s)
          b(c2) = b(c2) + A21*phi_face(s)
        end if
      else if(c2.lt.0) then

        ! Outflow is included because of the flux 
        ! corrections which also affects velocities
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) .or.  &
            (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)   .or.  &
            (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then    
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1)  = b(c1)  + A12 * phi % n(c2)

        ! Buffer: System matrix and parts belonging 
        ! to other subdomains are filled here.
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          if(grid % material(c1) == grid % material(c2)) then
            A % bou(c2) = -A12  ! cool parallel stuff
          else
            b(c1) = b(c1) + A12*phi_face(s)
          end if

        ! In case of wallflux 
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          Stot  = sqrt(  grid % sx(s)*grid % sx(s)  &
                       + grid % sy(s)*grid % sy(s)  &
                       + grid % sz(s)*grid % sz(s))
          b(c1) = b(c1) + Stot * phi % q(c2)
        endif 
      end if

    end if

  end do  ! through sides

  !-----------------------------!
  !                             !
  !   Temporal discretization   !
  !                             !
  !-----------------------------!
   
  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion == ADAMS_BASHFORTH) then 
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % d_o(c) - 0.5*phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion == CRANK_NICOLSON) then 
    do c = 1, grid % n_cells
      b(c)  = b(c) + 0.5*phi % d_o(c)
    end do  
  end if

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff == ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % c_o(c) - 0.5*phi % c_oo(c)
    end do 
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff == CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      if( (phi % c(c)+phi % c_o(c))  >= 0) then
        b(c)  = b(c) + 0.5*(phi % c(c) + phi % c_o(c))
      else
        A % val(A % dia(c)) = A % val(A % dia(c)) &
             - 0.5 * (phi % c(c) + phi % c_o(c)) / (phi % n(c)+1.e-6)
      end if
    end do
  end if
 
  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff == FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      if(phi % c(c) >= 0) then
        b(c)  = b(c) + phi % c(c)
      else
        A % val(A % dia(c)) = A % val(A % dia(c))  &
                            - phi % c(c)/(phi % n(c)+1.e-6)
      end if
    end do
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; Linear interpolation
  if(td_inertia == LINEAR) then
    do c = 1, grid % n_cells
      A0 = capacity * density * grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c)  = b(c) + A0*phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia == PARABOLIC) then
    do c = 1, grid % n_cells
      A0 = capacity * density * grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * A0
      b(c)  = b(c) + 2.0 * A0 * phi % o(c) - 0.5 * A0 * phi % oo(c)
    end do
  end if

  if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
     turbulence_model == HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant /= HYBRID) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = -0.22*Tsc(c) *&
                   (uu%n(c)*phi_x(c)+uv%n(c)*phi_y(c)+uw%n(c)*phi_z(c))
        u2uj_phij(c) = -0.22*Tsc(c)*&
                   (uv%n(c)*phi_x(c)+vv%n(c)*phi_y(c)+vw%n(c)*phi_z(c))
        u3uj_phij(c) = -0.22*Tsc(c)*&
                   (uw%n(c)*phi_x(c)+vw%n(c)*phi_y(c)+ww%n(c)*phi_z(c))
      end do
      call Grad_Mod_For_Phi(grid, u1uj_phij, 1, u1uj_phij_x, .true.)
      call Grad_Mod_For_Phi(grid, u2uj_phij, 2, u2uj_phij_y, .true.)
      call Grad_Mod_For_Phi(grid, u3uj_phij, 3, u3uj_phij_z, .true.)
      do c = 1, grid % n_cells
        b(c) = b(c) - (  u1uj_phij_x(c)  &
                       + u2uj_phij_y(c)  &
                       + u3uj_phij_z(c) ) * grid % vol(c)
      end do

      !------------------------------------------------------------------!
      !   Here we clean up transport equation from the false diffusion   !
      !------------------------------------------------------------------!
      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        Prt1 = 1.0/( 0.5882 + 0.228*(vis_t(c1)/(viscosity+1.0e-12)) - 0.0441*                  &
              (vis_t(c1)/(viscosity+1.0e-12))**2.0*(1.0 - exp(-5.165*( viscosity/(vis_t(c1)+1.0e-12) ))) )
        Prt2 = 1.0/( 0.5882 + 0.228*(vis_t(c2)/(viscosity+1.0e-12)) - 0.0441*                  &
              (vis_t(c2)/(viscosity+1.0e-12))**2.0*(1.0 - exp(-5.165*( viscosity/(vis_t(c2)+1.0e-12) ))) )

        Prt = fw(s)*Prt1 + (1.0-fw(s))*Prt2
        if(c2 > 0 .or.  &
           c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then
          phixS1 = fw(s)*phi_x(c1) + (1.0-fw(s))*phi_x(c2) 
          phiyS1 = fw(s)*phi_y(c1) + (1.0-fw(s))*phi_y(c2)
          phizS1 = fw(s)*phi_z(c1) + (1.0-fw(s))*phi_z(c2)
          phixS2 = phixS1 
          phiyS2 = phiyS1 
          phizS2 = phizS1 
          CONeff1 =       grid % f(s)  * (capacity*vis_t(c1)/Prt )  &
                  + (1. - grid % f(s)) * (capacity*vis_t(c2)/Prt )
          CONeff2 = CONeff1 
        else
          phixS1 = phi_x(c1) 
          phiyS1 = phi_y(c1) 
          phizS1 = phi_z(c1) 
          phixS2 = phixS1 
          phiyS2 = phiyS1 
          phizS2 = phizS1 
          CONeff1 = capacity*vis_t(c1)/Prt   
          CONeff2 = CONeff1 
        endif

        ! Total (exact) diffusive flux
        FUex1 = CONeff1 * (  phixS1*grid % sx(s)  &
                           + phiyS1*grid % sy(s)  &
                           + phizS1*grid % sz(s))
        FUex2 = CONeff2 * (  phixS2*grid % sx(s)  &
                           + phiyS2*grid % sy(s)  &
                           + phizS2*grid % sz(s))

        ! Implicit diffusive flux
        FUim1 = CONeff1*f_coef(s)*    &
                (  phixS1*grid % dx(s)      &
                 + phiyS1*grid % dy(s)      &
                 + phizS1*grid % dz(s) )
        FUim2 = CONeff2*f_coef(s)*    &
                (  phixS2*grid % dx(s)      &
                 + phiyS2*grid % dy(s)      &
                 + phizS2*grid % dz(s) )

        b(c1) = b(c1) - CONeff1*(phi%n(c2)-phi%n(c1))*f_coef(s)- FUex1 + FUim1
        if(c2  > 0) then
         b(c2) = b(c2) + CONeff1*(phi%n(c2)-phi%n(c1))*f_coef(s)+ FUex2 - FUim2
        end if
      end do
    end if  
  end if  

  call User_Mod_Source(grid, phi, A, b)

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
    call Control_Mod_Simple_Underrelaxation_For_Energy(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - urf) * phi % n(c)  &
                                      / urf
    A % val(A % dia(c)) = A % val(A % dia(c)) / urf
  end do  

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Energy_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of iterations based on coupling scheme
  if(coupling == 'PROJECTION') niter = 10
  if(coupling == 'SIMPLE')     niter =  5

  call bicg(A, phi % n, b, precond, niter, tol, ini_res, phi % res)

  call Info_Mod_Iter_Fill_At(2, 4, phi % name, niter, phi % res)

  call Comm_Mod_Exchange(grid, phi % n)

  end subroutine
