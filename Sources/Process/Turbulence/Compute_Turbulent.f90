!==============================================================================!
  subroutine Compute_Turbulent(grid, var, phi, Nstep)
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
  use Grid_Mod
  use Var_Mod
  use Work_Mod, only: phi_x => r_cell_01,  &
                      phi_y => r_cell_02,  &
                      phi_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: var
  type(Var_Type)  :: phi
  integer         :: Nstep
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, niter, miter, mat
  real    :: Fex, Fim 
  real    :: phis
  real    :: A0, A12, A21
  real    :: error
  real    :: VISeff
  real    :: phi_x_f, phi_y_f, phi_z_f
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
  call GraPhi(grid, phi % n, 1, phi_x, .TRUE.)
  call GraPhi(grid, phi % n, 2, phi_y, .TRUE.)
  call GraPhi(grid, phi % n, 3, phi_z, .TRUE.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Compute phimax and phimin
  do mat = 1, grid % n_materials
    if(BLEND_TUR(mat) /= NO) then
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
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      phis=f(s)*phi % n(c1) + (1.0-f(s))*phi % n(c2)

      if(BLEND_TUR(material(c1)) /= NO .or.  &
         BLEND_TUR(material(c2)) /= NO) then
        call Advection_Scheme(grid, phis, s, phi % n, phi_x, phi_y, phi_z,  &
                        grid % dx, grid % dy, grid % dz,        &
                        max(BLEND_TUR(material(c1)),            &
                            BLEND_TUR(material(c2))) ) 
      end if 

      ! Central differencing for advection
      if(ini == 1) then 
        if(c2  > 0) then
          phi % a_o(c1)=phi % a_o(c1)-Flux(s)*phis
          phi % a_o(c2)=phi % a_o(c2)+Flux(s)*phis
        else
          phi % a_o(c1)=phi % a_o(c1)-Flux(s)*phis
        endif 
      end if
      if(c2  > 0) then
        phi % a(c1)=phi % a(c1)-Flux(s)*phis
        phi % a(c2)=phi % a(c2)+Flux(s)*phis
      else
        phi % a(c1)=phi % a(c1)-Flux(s)*phis
      endif 

      ! Upwind 
      if(BLEND_TUR(material(c1)) /= NO .or. BLEND_TUR(material(c2)) /= NO) then
        if(Flux(s)  < 0) then   ! from c2 to c1
        phi % c(c1)=phi % c(c1)-Flux(s)*phi % n(c2)
          if(c2  > 0) then
            phi % c(c2)=phi % c(c2)+Flux(s)*phi % n(c2)
          endif
        else 
          phi % c(c1)=phi % c(c1)-Flux(s)*phi % n(c1)
          if(c2  > 0) then
            phi % c(c2)=phi % c(c2)+Flux(s)*phi % n(c1)
          endif
        end if
      end if   ! BLEND_TUR 
    else       ! c2 < 0

    !=================================================
    !   I excluded this for the inflow-outflow 
    !-------------------------------------------------
    !   CUi(c1) = CUi(c1)-Flux(s)*phi(c2) corrected !
    !---- Full upwind
    !      if(TypeBC(c2) == OUTFLOW) then
    !      if(BLEND_TUR == YES) then
    !        if(Flux(s)  < 0) then   ! from c2 to c1
    !          XUi(c1)=XUi(c1)-Flux(s)*phi(c2)
    !        else 
    !          XUi(c1)=XUi(c1)-Flux(s)*phi(c1)
    !        end if 
    !      end if ! BLEND_TUR = YES
    !      end if   ! OUTFLOW
    !=================================================
    end if     ! c2 > 0 
  end do    ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(CONVEC == AB) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (1.5*phi % a_o(c) - 0.5*phi % a_oo(c) - phi % c(c))
    end do  
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(CONVEC == CN) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (0.5 * ( phi % a(c) + phi % a_o(c) ) - phi % c(c))
    end do  
  endif

  ! Fully implicit treatment of convective fluxes 
  if(CONVEC == FI) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (phi % a(c) - phi % c(c))
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

    VISeff = VISc + (fF(s)*VISt(c1) + (1.0-fF(s))*VISt(c2))/phi % Sigma 

    if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA)          &
    VISeff = VISc+(fF(s)*VIS % n(c1)+(1.0-fF(s))*VIS % n(c2))/phi % Sigma

    if(SIMULA==HYB_ZETA)          &
    VISeff = VISc + (fF(s)*VISt_eff(c1) + (1.0-fF(s))*VISt_eff(c2))/phi % Sigma

    phi_x_f = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2)
    phi_y_f = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
    phi_z_f = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)

    ! This implements zero gradient for k
    if(SIMULA==K_EPS.and.MODE==HIGH_RE) then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then  
          phi_x_f = 0.0 
          phi_y_f = 0.0
          phi_z_f = 0.0
          VISeff  = 0.0
        end if 
      end if
    end if

    if(SIMULA==ZETA.or.SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA) then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
          if(sqrt(TauWall(c1))*WallDs(c1)/VISc>2.0) then      
            phi_x_f = 0.0
            phi_y_f = 0.0
            phi_z_f = 0.0
            VISeff  = 0.0
          end if
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    Fex = VISeff * (  phi_x_f * grid % sx(s)  &
                    + phi_y_f * grid % sy(s)  &
                    + phi_z_f * grid % sz(s) )

    A0 = VISeff * Scoef(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: Scoef is
    !  not corrected at interface between materials)
    Fim = (  phi_x_f * grid % dx(s)                      &
           + phi_y_f * grid % dy(s)                      &
           + phi_z_f * grid % dz(s) ) * A0

    ! This is yet another crude approximation:
    ! A0 is calculated approximatelly
    !    if( StateMat(material(c1))==FLUID .and.  &  ! 2mat
    !        StateMat(material(c2))==SOLID        &  ! 2mat
    !        .or.                                 &  ! 2mat 
    !        StateMat(material(c1))==SOLID .and.  &  ! 2mat
    !        StateMat(material(c2))==FLUID ) then    ! 2mat
    !      A0 = A0 + A0                              ! 2mat
    !    end if                                      ! 2mat

    ! Straight diffusion part 
    if(ini == 1) then
      if(c2  > 0) then
        phi % d_o(c1) = phi % d_o(c1) + (phi % n(c2)-phi % n(c1))*A0   
        phi % d_o(c2) = phi % d_o(c2) - (phi % n(c2)-phi % n(c1))*A0    
      else
        if(TypeBC(c2) /= SYMMETRY) then
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
    if( (DIFFUS == CN) .or. (DIFFUS == FI) ) then  
      if(DIFFUS  ==  CN) then       ! Crank Nicholson
        A12 = 0.5 * A0 
        A21 = 0.5 * A0 
      end if

      if(DIFFUS  ==  FI) then       ! Fully implicit
        A12 = A0 
        A21 = A0
      end if

      if(BLEND_TUR(material(c1)) /= NO .or. BLEND_TUR(material(c2)) /= NO) then
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
        if((TypeBC(c2) == INFLOW)  .or.                 &
           (TypeBC(c2) == WALL)    .or.                 &
           (TypeBC(c2) == PRESSURE).or.                 &
           (TypeBC(c2) == CONVECT) .or.                 &
           (TypeBC(c2) == WALLFL) ) then                               
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1) = b(c1) + A12 * phi % n(c2)
        else if( TypeBC(c2) == BUFFER ) then  
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
  if(DIFFUS == AB) then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % d_o(c) - 0.5 * phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % d_o(c)
    end do  
  end if
             
  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % c_o(c) - 0.5 * phi % c_oo(c)
    end do 
  end if
    
  ! Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % c(c) + 0.5 * phi % c_o(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
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
  if(INERT == LIN) then
    do c = 1, grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(INERT == PAR) then
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
  if(SIMULA == K_EPS) then 
    if(phi % name == 'KIN') call Source_Kin_K_Eps(grid)
    if(phi % name == 'EPS') call Source_Eps_K_Eps(grid)
  end if

  if(SIMULA == K_EPS_VV .or.  &
     SIMULA == ZETA     .or.  &
     SIMULA == HYB_ZETA) then
    if(phi % name == 'KIN') call Source_Kin_K_Eps_V2_F(grid)
    if(phi % name == 'EPS') call Source_Eps_K_Eps_V2_F(grid)
    if(phi % name == 'V^2') call Source_V2_K_Eps_V2_F(grid, Nstep)
  end if

  if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA)                                &
  call Source_Vis_Spalart_Almaras(grid, phi_x, phi_y, phi_z)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + A % val(A % dia(c)) * (1.0-phi % URF)*phi % n(c) / phi % URF
    A % val(A % dia(c)) = A % val(A % dia(c)) / phi % URF
  end do

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5

  niter=miter
  call cg(A, phi % n, b,            &
          PREC, niter, phi % STol,  &
          res(var), error)

  do c = 1, grid % n_cells
    if( phi%n(c)<0.0)then
      phi%n(c) = phi%o(c)
    end if
  end do 

  if(this_proc < 2) write(*,*) '# ', phi % name, res(var), niter 

  call Exchange(grid, phi % n)

  end subroutine
