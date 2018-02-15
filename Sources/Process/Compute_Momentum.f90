!==============================================================================!
  subroutine Compute_Momentum(grid, dt, ini, var, ui,                      &
                              ui_i, ui_j, ui_k,                            &
                              Si, Sj, Sk,                                  &
                              Di, Dj, Dk,                                  &
                              Hi, uj_i, uk_i)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
  use Var_Mod
  use Grid_Mod
  use Bulk_Mod
  use Info_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: ini
  real            :: dt
  integer         :: var
  type(Var_Type)  :: ui
  real            :: ui_i(-grid % n_bnd_cells:grid % n_cells),  &
                     ui_j(-grid % n_bnd_cells:grid % n_cells),  & 
                     ui_k(-grid % n_bnd_cells:grid % n_cells)
  real            :: Si(grid % n_faces),  &
                     Sj(grid % n_faces),  &
                     Sk(grid % n_faces) 
  real            :: Di(grid % n_faces),  &
                     Dj(grid % n_faces),  &
                     Dk(grid % n_faces) 
  real            :: Hi  (-grid % n_bnd_cells:grid % n_cells),  &
                     uj_i(-grid % n_bnd_cells:grid % n_cells),  &
                     uk_i(-grid % n_bnd_cells:grid % n_cells) 
  real            :: uuS, vvS, wwS, uvS, uwS, vwS
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter, miter, mat
  real              :: Fex, Fim 
  real              :: uis
  real              :: A0, A12, A21
  real              :: error, tol
  real              :: VISeff, vis_tS, Fstress 
  real              :: ui_iS,ui_jS,ui_kS,uj_iS,uk_iS
  character(len=80) :: coupling
  character(len=80) :: precond
  character(len=80) :: sd_advection  ! space disretization of advection (scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  character(len=80) :: td_inertia    ! time-disretization for inerita  
  character(len=80) :: td_advection  ! time-disretization for advection
  character(len=80) :: td_diffusion  ! time-disretization for diffusion 
  character(len=80) :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor                 
  character(len=80) :: heat_transfer
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
!------------------------------------------------------------------------------!
!
!  Stress tensor on the face s:
!
!    T = mi * [    2*dU/dx     dU/dy+dV/dx   dU/dz+dW/dx  ]
!             [  dU/dy+dV/dx     2*dV/dy     dV/dz+dW/dy  ]
!             [  dU/dz+dW/dx   dV/dz+dW/dy     2*dW/dz    ]
!
!  The forces, acting on the cell face are:
!  
!    Fx = T11*Sx + T12*Sy + T13*Sz
!    Fy = T21*Sx + T22*Sy + T23*Sz
!    Fz = T31*Sx + T32*Sy + T33*Sz
!
!  which could also be written in the compact form:
!
!    {F} = [T]{S}
!
!  U component:
!
!    Fx = Txx*Sx + Txy*Sy + Txz*Sz 
!
!  V component:    
!
!    Fy = Tyx*Sx + Tyy*Sy + Tyz*Sz   
!
!  W component:
!
!    Fz = Tzx*Sx + Tzy*Sy + Tzz*Sz
!
!------------------------------------------------------------------------------!
!     
!  The form of equations which I am solving:    
!     
!     /             /              /               /
!    |     du      |              |               |
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | p dS
!    |     dt      |              |               |
!   /             /              /               /
!
!  Dimension of the system under consideration
!
!     [A]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     A              [kg/s]
!     U, V, W        [m/s]
!     bU, bV, bW     [kgm/s^2]   [N]
!     P, PP          [kg/ms^2]   [N/m^2]
!     Flux           [kg/s]
!     CU*, CV*, CW*  [kgm/s^2]   [N]
!     DU*, DV*, DW*  [kgm/s^2]   [N]
!     XU*, XV*, XW*  [kgm/s^2]   [N]
!
!==============================================================================!

  b = 0.0
  A % val = 0.0
  Fstress = 0.0

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

  call Control_Mod_Heat_Transfer(heat_transfer)
  call Control_Mod_Turbulence_Model(turbulence_model)
  call Control_Mod_Turbulence_Model_Variant(turbulence_model_variant)

  ! Old values (o) and older than old (oo)
  if(ini == 1) then
    do c=1,grid % n_cells
      ui % oo(c)  = ui % o(c)
      ui % o (c)  = ui % n(c)
      ui % a_oo(c) = ui % a_o(c)
      ui % a_o (c) = 0.0 
      ui % d_oo(c) = ui % d_o(c)
      ui % d_o (c) = 0.0 
      ui % c_oo(c) = ui % c_o(c)
      ui % c_o (c) = ui % c(c) 
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Momentum(sd_advection)
  call Control_Mod_Blending_Coefficient_Momentum(blend)
  
  ! Compute phimax and phimin
  do mat=1,grid % n_materials
    if(sd_advection .ne. 'CENTRAL') then
      call Compute_Minimum_Maximum(grid, ui % n)  ! or ui % o ???
      goto 1  ! why on Earth this?
    end if
  end do

  ! New values
1 do c=1,grid % n_cells
    ui % a(c)    = 0.0
    ui % c(c)    = 0.0
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s=1,grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s) 

    ! Central differencing
    uis = grid % f(s) * ui % n(c1) + (1.0 - grid % f(s)) * ui % n(c2)

    if(sd_advection .ne. 'CENTRAL') then
      call Advection_Scheme(grid, uis, s, ui % n,    &
                            ui_i, ui_j, ui_k,        &
                            Di, Dj, Dk,              &
                            sd_advection, blend) 
    end if 
    
    ! Compute advection term
    if(ini == 1) then 
      if(c2  > 0) then
        ui % a_o(c1)=ui % a_o(c1)-Flux(s)*uis
        ui % a_o(c2)=ui % a_o(c2)+Flux(s)*uis
      else
        ui % a_o(c1)=ui % a_o(c1)-Flux(s)*uis
      endif 
    end if

    if(c2  > 0) then
      ui % a(c1)=ui % a(c1)-Flux(s)*uis
      ui % a(c2)=ui % a(c2)+Flux(s)*uis
    else
      ui % a(c1)=ui % a(c1)-Flux(s)*uis
    endif 

    ! Store upwinded part of the advection term in "c"
    if(coupling .ne. 'PROJECTION') then
      if(Flux(s)  < 0) then   ! from c2 to c1
        ui % c(c1)=ui % c(c1)-Flux(s)*ui % n(c2)
        if(c2  > 0) then
          ui % c(c2)=ui % c(c2)+Flux(s)*ui % n(c2)
        endif
      else 
        ui % c(c1)=ui % c(c1)-Flux(s)*ui % n(c1)
        if(c2  > 0) then
          ui % c(c2)=ui % c(c2)+Flux(s)*ui % n(c1)
        endif
      end if
    end if   ! BLEND 
  end do    ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(td_advection == 'ADAMS_BASHFORTH') then
    do c=1,grid % n_cells
      b(c) = b(c) + (1.5*ui % a_o(c) - 0.5*ui % a_oo(c) - ui % c(c))
    end do  
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(td_advection == 'CRANK_NICOLSON') then
    do c=1,grid % n_cells
      b(c) = b(c) + (0.5 * ( ui % a(c) + ui % a_o(c) ) - ui % c(c))
    end do  
  endif

  ! Fully implicit treatment of convective fluxes 
  if(td_advection == 'FULLY_IMPLICIT') then
    do c=1,grid % n_cells
      b(c) = b(c) + (ui % a(c) - ui % c(c))
    end do  
  end if     
          
  ! New values
  do c=1,grid % n_cells
    ui % c(c) = 0.0
  end do

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s=1,grid % n_faces       

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)   

    VISeff = fF(s)*vis_t(c1)+(1.0-fF(s))*vis_t(c2) + VISc

    if(turbulence_model=='HYB_ZETA') then
      VISeff = fF(s)*vis_t_eff(c1)+(1.0-fF(s))*vis_t_eff(c2) + VISc
    end if

    if(c2 < 0 .and. turbulence_model == 'LES') then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
        VISeff = VISwall(c1)
      end if
    end if 

    if( turbulence_model == 'ZETA'             .or.  &
        turbulence_model == 'K_EPS_VV'         .or.  &
       (turbulence_model == 'K_EPS' .and.  &
        turbulence_model_variant == 'HIGH_RE') .or.  &
        turbulence_model == 'HYB_ZETA') then
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then
          VISeff = VISwall(c1)
        end if
      end if
    end if

    ! Add influence of Re stresses for 'EBM'
    if(turbulence_model == 'EBM'.or.turbulence_model == 'HJ') then
      if(turbulence_model_variant /= 'HYBRID') then        
        if(ui % name == 'U') then
          uuS = fF(s)*uu % n(c1)+(1.0-fF(s))*uu % n(c2)
          uvS = fF(s)*uv % n(c1)+(1.0-fF(s))*uv % n(c2)
          uwS = fF(s)*uw % n(c1)+(1.0-fF(s))*uw % n(c2)
          Fstress = - (  uuS * grid % sx(s)  &
                       + uvS * grid % sy(s)  &
                       + uwS * grid % sz(s) )  
        else if(ui % name == 'V') then 
          uvS = fF(s)*uv % n(c1)+(1.0-fF(s))*uv % n(c2)
          vvS = fF(s)*vv % n(c1)+(1.0-fF(s))*vv % n(c2)
          vwS = fF(s)*vw % n(c1)+(1.0-fF(s))*vw % n(c2)
          Fstress =  - (  uvS * grid % sx(s)  &
                        + vvS * grid % sy(s)  &
                        + vwS * grid % sz(s) )  
        else if(ui % name == 'W') then 
          uwS = fF(s)*uw % n(c1)+(1.0-fF(s))*uw % n(c2)
          vwS = fF(s)*vw % n(c1)+(1.0-fF(s))*vw % n(c2)
          wwS = fF(s)*ww % n(c1)+(1.0-fF(s))*ww % n(c2)
          Fstress =  - (  uwS * grid % sx(s)  &
                        + vwS * grid % sy(s)  &
                        + wwS * grid % sz(s) )  
        end if 
      end if 
    end if

    ui_iS = fF(s)*ui_i(c1) + (1.0-fF(s))*ui_i(c2)
    ui_jS = fF(s)*ui_j(c1) + (1.0-fF(s))*ui_j(c2)
    ui_kS = fF(s)*ui_k(c1) + (1.0-fF(s))*ui_k(c2)
    uj_iS = fF(s)*uj_i(c1) + (1.0-fF(s))*uj_i(c2)
    uk_iS = fF(s)*uk_i(c1) + (1.0-fF(s))*uk_i(c2)

    ! total (exact) viscous stress 
    Fex=VISeff*(      2.0*ui_iS   * Si(s)      &
                  + (ui_jS+uj_iS) * Sj(s)      & 
                  + (ui_kS+uk_iS) * Sk(s) )

    A0 = VISeff * Scoef(s)

    ! Implicit viscous stress
    ! this is a very crude approximation: Scoef is not
    ! corrected at interface between materials
    Fim=(   ui_iS*Di(s)                &
          + ui_jS*Dj(s)                &
          + ui_kS*Dk(s))*A0

    ! Straight diffusion part 
    if(ini == 1) then
      if(c2  > 0) then
        ui % d_o(c1) = ui % d_o(c1) + (ui % n(c2)-ui % n(c1))*A0   
        ui % d_o(c2) = ui % d_o(c2) - (ui % n(c2)-ui % n(c1))*A0    
      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) /= SYMMETRY) then
          ui % d_o(c1) = ui % d_o(c1) + (ui % n(c2)-ui % n(c1))*A0   
        end if 
      end if 
    end if

    ! Cross diffusion part
    ui % c(c1) = ui % c(c1) + Fex - Fim + Fstress
    if(c2  > 0) then
      ui % c(c2) = ui % c(c2) - Fex + Fim - Fstress
    end if 

     ! Compute the coefficients for the sysytem matrix
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
        if((Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW)  .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL)    .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT) .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL)) then                                
           ! (Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW) ) then   
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1) = b(c1) + A12 * ui % n(c2)
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then  
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          A % bou(c2) = -A12  ! cool parallel stuff
        endif
      end if     
    end if
  end do  ! through faces

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  !---------------------------------!
  !
  ! Add Re stress influence on momentum
  ! This is an alternative way to implement RSM. 
  !
  !  if(turbulence_model == 'EBM'.or.turbulence_model == 'HJ') then
  !    if(ui % name == 'U') then
  !      call GraPhi(uu%n,1,VAR2x,.TRUE.)
  !      call GraPhi(uv%n,2,VAR2y,.TRUE.)
  !      call GraPhi(uw%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    else if(ui % name == 'V') then
  !      call GraPhi(uv%n,1,VAR2x,.TRUE.)
  !      call GraPhi(vv%n,2,VAR2y,.TRUE.)
  !      call GraPhi(vw%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    else if(ui % name == 'W') then
  !      call GraPhi(uw%n,1,VAR2x,.TRUE.)
  !      call GraPhi(vw%n,2,VAR2y,.TRUE.)
  !      call GraPhi(ww%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    end if
  
  ! Here we clean up momentum from the false diffusion
  if(turbulence_model == 'EBM' .or.  &
     turbulence_model == 'HJ') then
    if(turbulence_model_variant /= 'HYBRID') then        
      do s=1,grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        vis_tS = (fF(s)*vis_t(c1)+(1.0-fF(s))*vis_t(c2))
        A0 = Scoef(s)*vis_tS 
        VISeff = vis_tS

        ui_iS = fF(s)*ui_i(c1) + (1.0-fF(s))*ui_i(c2)
        ui_jS = fF(s)*ui_j(c1) + (1.0-fF(s))*ui_j(c2)
        ui_kS = fF(s)*ui_k(c1) + (1.0-fF(s))*ui_k(c2)
        uj_iS = fF(s)*uj_i(c1) + (1.0-fF(s))*uj_i(c2)
        uk_iS = fF(s)*uk_i(c1) + (1.0-fF(s))*uk_i(c2)

        Fex=VISeff*( 2.0*ui_iS*Si(s)                           &
                      + (ui_jS+uj_iS)*Sj(s)                   &
                      + (ui_kS+uk_iS)*Sk(s) )

        Fim=( ui_iS*Di(s)                                      &
             +ui_jS*Dj(s)                                      &
             +ui_kS*Dk(s))*VISeff*Scoef(s)

        b(c1) = b(c1) - VISeff*(ui%n(c2)-ui%n(c1))*Scoef(s)- Fex + Fim
        if(c2  > 0) then
          b(c2) = b(c2) + VISeff*(ui%n(c2)-ui%n(c1))*Scoef(s)+ Fex - Fim
        end if
      end do
    end if 
  end if

  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion == 'ADAMS_BASHFORTH') then 
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * ui % d_o(c) - 0.5 * ui % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion == 'CRANK_NICOLSON') then 
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * ui % d_o(c)
    end do  
  end if
             
  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff == 'ADAMS_BASHFORTH') then
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * ui % c_o(c) - 0.5 * ui % c_oo(c)
    end do 
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff == 'CRANK_NICOLSON') then
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * ui % c(c) + 0.5 * ui % c_o(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff == 'FULLY_IMPLICIT') then
    do c=1,grid % n_cells
      b(c) = b(c) + ui % c(c)
    end do 
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; linear interpolation
  if(td_inertia == 'LINEAR') then
    do c=1,grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * ui % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia == 'PARABOLIC') then
    do c=1,grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * A0
      b(c) = b(c) + 2.0*A0 * ui % o(c) - 0.5*A0 * ui % oo(c)
    end do
  end if

  !---------------------------------!
  !                                 !
  !   Pressure term contributions   !
  !                                 !
  !---------------------------------!

  !--------------------------!
  !   Global pressure drop   !
  !--------------------------!
  if(ui % name == 'U') then
    do c=1,grid % n_cells
      b(c) = b(c) + bulk(material(c)) % p_drop_x * grid % vol(c)
    end do
  else if(ui % name == 'V') then
    do c=1,grid % n_cells
      b(c) = b(c) + bulk(material(c)) % p_drop_y * grid % vol(c)
    end do
  else if(ui % name == 'W') then
    do c=1,grid % n_cells
      b(c) = b(c) + bulk(material(c)) % p_drop_z * grid % vol(c)
    end do
  end if

  !---------------------------------!
  !   Local pressure distribution   !
  !---------------------------------!
  do c=1,grid % n_cells
    b(c) = b(c) - Hi(c) * grid % vol(c)
  end do

  !----------------------------------------!
  !   All other terms defined by the user  !
  !----------------------------------------!
  if(heat_transfer == 'YES') call User_Mod_Force(grid, ui, a, b)

  !-----------------------------------!
  !                                   !
  !   Solve the equations for u,v,w   !
  !                                   !    
  !-----------------------------------!

  ! Type of coupling is important
  call Control_Mod_Pressure_Momentum_Coupling(coupling)

  ! Set under-relaxation factor
  urf = 1.0
  if(coupling == 'SIMPLE')  &
    call Control_Mod_Simple_Underrelaxation_For_Momentum(urf)

  do c=1,grid % n_cells
    A % sav(c) = A % val(A % dia(c))
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - urf)*ui % n(c) / urf
    A % val(A % dia(c)) = A % val(A % dia(c)) / urf
  end do

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Momentum_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of solver iterations on coupling method
  if(coupling == 'PROJECTION') miter = 10
  if(coupling == 'SIMPLE')     miter =  5

  niter=miter

  call cg(A, ui % n, b, precond, niter, tol, res(var), error)

  if(ui % name == 'U') then
    call Info_Mod_Iter_Fill_At(2, 1, ui % name, niter, res(var))
  end if
  if(ui % name == 'V') then
    call Info_Mod_Iter_Fill_At(2, 2, ui % name, niter, res(var))
  end if
  if(ui % name == 'W') then
    call Info_Mod_Iter_Fill_At(2, 3, ui % name, niter, res(var))
  end if

  call Exchange(grid, ui % n)

  end subroutine
