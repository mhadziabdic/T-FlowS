!==============================================================================!
  subroutine Compute_Turbulent(grid, var, phi, phi_x, phi_y, phi_z, Nstep)
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
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: var
  type(Var_Type)  :: phi
  real            :: phi_x(-NbC:NC), phi_y(-NbC:NC), phi_z(-NbC:NC)
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, niter, miter, mat, Nstep
  real    :: Fex, Fim 
  real    :: phis
  real    :: A0, A12, A21
  real    :: error
  real    :: VISeff
  real    :: phi_xS, phi_yS, phi_zS
!==============================================================================!
!                                                                              ! 
!  The form of equations which are solved:                                     !   
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!------------------------------------------------------------------------------!

  A % val = 0.0

  b=0.0

  ! This is important for "copy" boundary conditions. Find out why !
  do c=-NbC,-1
    A % bou(c)=0.0
  end do

  !-------------------------------------! 
  !   Initialize variables and fluxes   !
  !-------------------------------------! 

  ! Old values (o) and older than old (oo)
  if(ini == 1) then
    do c=1,NC
      phi % oo(c)  = phi % o(c)
      phi % o (c)  = phi % n(c)
      phi % Coo(c) = phi % Co(c)
      phi % Co (c) = 0.0 
      phi % Doo(c) = phi % Do(c)
      phi % Do (c) = 0.0 
      phi % Xoo(c) = phi % Xo(c)
      phi % Xo (c) = phi % X(c) 
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Compute phimax and phimin
  do mat = 1, grid % n_materials
    if(BLEND_TUR(mat) /= NO) then
      call CalMinMax(phi % n)  ! or phi % o ???
      goto 1
    end if
  end do

  ! New values
1 do c=1,NC
    phi % C(c)    = 0.0
    phi % X(c)    = 0.0
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s=1,NS

    c1=SideC(1,s)
    c2=SideC(2,s) 

    ! Velocities on "orthogonal" cell centers 
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      phis=f(s)*phi % n(c1) + (1.0-f(s))*phi % n(c2)

      if(BLEND_TUR(material(c1)) /= NO .or.  &
         BLEND_TUR(material(c2)) /= NO) then
        call ConvScheme(phis, s, phi % n, phi_x, phi_y, phi_z,  &
                        grid % dx, grid % dy, grid % dz,        &
                        max(BLEND_TUR(material(c1)),            &
                            BLEND_TUR(material(c2))) ) 
      end if 

      ! Central differencing for advection
      if(ini == 1) then 
        if(c2  > 0) then
          phi % Co(c1)=phi % Co(c1)-Flux(s)*phis
          phi % Co(c2)=phi % Co(c2)+Flux(s)*phis
        else
          phi % Co(c1)=phi % Co(c1)-Flux(s)*phis
        endif 
      end if
      if(c2  > 0) then
        phi % C(c1)=phi % C(c1)-Flux(s)*phis
        phi % C(c2)=phi % C(c2)+Flux(s)*phis
      else
        phi % C(c1)=phi % C(c1)-Flux(s)*phis
      endif 

      ! Upwind 
      if(BLEND_TUR(material(c1)) /= NO .or. BLEND_TUR(material(c2)) /= NO) then
        if(Flux(s)  < 0) then   ! from c2 to c1
        phi % X(c1)=phi % X(c1)-Flux(s)*phi % n(c2)
          if(c2  > 0) then
            phi % X(c2)=phi % X(c2)+Flux(s)*phi % n(c2)
          endif
        else 
          phi % X(c1)=phi % X(c1)-Flux(s)*phi % n(c1)
          if(c2  > 0) then
            phi % X(c2)=phi % X(c2)+Flux(s)*phi % n(c1)
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
    do c=1,NC
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (1.5*phi % Co(c) - 0.5*phi % Coo(c) - phi % X(c))
    end do  
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(CONVEC == CN) then
    do c=1,NC
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (0.5 * ( phi % C(c) + phi % Co(c) ) - phi % X(c))
    end do  
  endif

  ! Fully implicit treatment of convective fluxes 
  if(CONVEC == FI) then
    do c=1,NC
      b(c) = b(c) + URFC_Tur(material(c)) * &
                    (phi % C(c) - phi % X(c))
    end do  
  end if     
          
  ! New values
  do c=1,NC
    phi % X(c) = 0.0
  end do
  
  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s=1,NS       

    c1=SideC(1,s)
    c2=SideC(2,s)   

    VISeff = VISc + (fF(s)*VISt(c1) + (1.0-fF(s))*VISt(c2))/phi % Sigma 

    if(SIMULA==SPA_ALL.or.SIMULA==DES_SPA)          &
    VISeff = VISc+(fF(s)*VIS % n(c1)+(1.0-fF(s))*VIS % n(c2))/phi % Sigma

    if(SIMULA==HYB_ZETA)          &
    VISeff = VISc + (fF(s)*VISt_eff(c1) + (1.0-fF(s))*VISt_eff(c2))/phi % Sigma

    phi_xS = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2)
    phi_yS = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
    phi_zS = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)

    ! This implements zero gradient for k
    if(SIMULA==K_EPS.and.MODE==HRe) then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then  
          phi_xS = 0.0 
          phi_yS = 0.0
          phi_zS = 0.0
          VISeff  = 0.0
        end if 
      end if
    end if

    if(SIMULA==ZETA.or.SIMULA==K_EPS_VV.or.SIMULA==HYB_ZETA) then
      if(c2 < 0 .and. phi % name == 'KIN') then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
          if(sqrt(TauWall(c1))*WallDs(c1)/VISc>2.0) then      
            phi_xS = 0.0
            phi_yS = 0.0
            phi_zS = 0.0
            VISeff  = 0.0
          end if
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    Fex = VISeff * (  phi_xS * grid % sx(s)  &
                    + phi_yS * grid % sy(s)  &
                    + phi_zS * grid % sz(s) )

    A0 = VISeff * Scoef(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: Scoef is
    !  not corrected at interface between materials)
    Fim = (  phi_xS * grid % dx(s)                      &
           + phi_yS * grid % dy(s)                      &
           + phi_zS * grid % dz(s) ) * A0

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
        phi % Do(c1) = phi % Do(c1) + (phi % n(c2)-phi % n(c1))*A0   
        phi % Do(c2) = phi % Do(c2) - (phi % n(c2)-phi % n(c1))*A0    
      else
        if(TypeBC(c2) /= SYMMETRY) then
          phi % Do(c1) = phi % Do(c1) + (phi % n(c2)-phi % n(c1))*A0   
        end if 
      end if 
    end if

    ! Cross diffusion part
    phi % X(c1) = phi % X(c1) + Fex - Fim 
    if(c2  > 0) then
      phi % X(c2) = phi % X(c2) - Fex + Fim 
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
    do c=1,NC
      b(c) = b(c) + 1.5 * phi % Do(c) - 0.5 * phi % Doo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c=1,NC
      b(c) = b(c) + 0.5 * phi % Do(c)
    end do  
  end if
             
  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c=1,NC
      b(c) = b(c) + 1.5 * phi % Xo(c) - 0.5 * phi % Xoo(c)
    end do 
  end if
    
  ! Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c=1,NC
      b(c) = b(c) + 0.5 * phi % X(c) + 0.5 * phi % Xo(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
    do c=1,NC
      b(c) = b(c) + phi % X(c)
    end do 
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; linear interpolation
  if(INERT == LIN) then
    do c=1,NC
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(INERT == PAR) then
    do c=1,NC
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
  do c=1,NC
    b(c) = b(c) + A % val(A % dia(c)) * (1.0-phi % URF)*phi % n(c) / phi % URF
    A % val(A % dia(c)) = A % val(A % dia(c)) / phi % URF
!?????? Asave(c) = A % val(A % dia(c)) ??????
  end do

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5

  niter=miter
  call cg(NC, Nbc, A,           &
         phi % n, b, PREC,      &
         niter,phi % STol, res(var), error)

  do c= 1, NC
    if( phi%n(c)<0.0)then
      phi%n(c) = phi%o(c)
    end if
  end do 

  if(this_proc < 2) write(*,*) '# ', phi % name, res(var), niter 

  call Exchng(phi % n)

  end subroutine
