!==============================================================================!
  subroutine Compute_Stresses(grid, var, phi,             &
                              phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!   EBM and HJ are calling this subroutine.                                    !
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
  integer :: s, c, c1, c2, niter, miter, mat
  real    :: Fex, Fim, TDiffS, VIStur 
  real    :: phis
  real    :: A0, A12, A21
  real    :: error
  real    :: VISeff
  real    :: phi_xS, phi_yS, phi_zS
  real    :: VIStS, TDiff_im, TDiff_ex, TDiff_coef, TDiffx, TDiffy,TDiffz
  real    :: uuS, vvS, wwS, uvS, uwS, vwS
!------------------------------------------------------------------------------!
!                                                                              ! 
!  The form of equations which are being solved:                               !   
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!==============================================================================!


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
  do mat=1,grid % n_materials
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

      if(BLEND_TUR(material(c1)) /= NO .or. BLEND_TUR(material(c2)) /= NO) then
        call ConvScheme(phis, s, phi % n, phi_x, phi_y, phi_z, Dx, Dy, Dz, &
                        max(BLEND_TUR(material(c1)),BLEND_TUR(material(c2))) ) 
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
    !      ! Full upwind
    !      if(TypeBC(c2) == OUTFLOW) then
    !        if(BLEND_TUR == YES) then
    !          if(Flux(s)  < 0) then   ! from c2 to c1
    !            XUi(c1)=XUi(c1)-Flux(s)*phi(c2)
    !          else 
    !            XUi(c1)=XUi(c1)-Flux(s)*phi(c1)
    !          end if 
    !        end if ! BLEND_TUR = YES
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

    ! VIStur is used to make diaginal element more dominant.
    ! This contribution is later substracted.
    VIStS = fF(s)*VISt(c1) + (1.0-fF(s))*VISt(c2)

    if(SIMULA==EBM) then
      VISeff = VISc + VIStS 
    else if(SIMULA==HJ.and.MODE/=HYB) then
      VISeff = 1.5*VISc 
    else if(SIMULA==HJ.and.MODE==HYB) then
      VISeff = VISc 
    end if

    if(SIMULA==HJ) then
      if(MODE==HYB) then
        VISeff = VISeff + VIStS
      end if
    end if

    phi_xS = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2)
    phi_yS = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
    phi_zS = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)


    ! Total (exact) diffusive flux plus turb. diffusion
    Fex=VISeff*(phi_xS*Sx(s) + phi_yS*Sy(s) + phi_zS*Sz(s)) 

    A0 = VISeff * Scoef(s) 

    ! Implicit diffusive flux
    ! (this is a very crude approximation: Scoef is
    !  not corrected at interface between materials)
    Fim=( phi_xS*Dx(s)                      &
         +phi_yS*Dy(s)                      &
         +phi_zS*Dz(s))*A0

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

      if(DIFFUS  ==  CN) then  ! Crank-Nicholson
        A12 = 0.5 * A0 
        A21 = 0.5 * A0 
      end if

      if(DIFFUS  ==  FI) then  ! fully implicit
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
        ! Convect is commented because for turbulent scalars convect 
        ! outflow is treated as classic outflow.
        if((TypeBC(c2) == INFLOW).or.                     &
           (TypeBC(c2) == WALL).or.                       &
!!!           (TypeBC(c2) == CONVECT).or.                    &
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

  !------------------------------!
  !   Turbulent diffusion term   !
  !------------------------------!
  if(MODE/=HYB) then
    if(var==13) then
      CmuD = 0.18        
    else
      CmuD = 0.22
    end if 
    if(SIMULA==HJ) then        
      do c=1,NC
        VAR1x(c) = CmuD/phi % Sigma*Kin%n(c)/Eps%n(c)*&
                   ((uu%n(c))*phi_x(c)+uv%n(c)*phi_y(c)+uw%n(c)*phi_z(c)) - &
                   VISc*phi_x(c)

        VAR1y(c) = CmuD/phi % Sigma*Kin%n(c)/Eps%n(c)*&
                   (uv%n(c)*phi_x(c)+(vv%n(c))*phi_y(c)+vw%n(c)*phi_z(c)) - &
                   VISc*phi_y(c)

        VAR1z(c) = CmuD/phi % Sigma*Kin%n(c)/Eps%n(c)*&
                   (uw%n(c)*phi_x(c)+vw%n(c)*phi_y(c)+(ww%n(c))*phi_z(c)) - &
                   VISc*phi_z(c) 
      end do
    else if(SIMULA==EBM) then
      do c=1,NC
        VAR1x(c) = CmuD/phi % Sigma*Tsc(c)*&
                   ((uu%n(c))*phi_x(c)+uv%n(c)*phi_y(c)+uw%n(c)*phi_z(c)) 

        VAR1y(c) = CmuD/phi % Sigma*Tsc(c)*&
                   (uv%n(c)*phi_x(c)+(vv%n(c))*phi_y(c)+vw%n(c)*phi_z(c)) 

        VAR1z(c) = CmuD/phi % Sigma*Tsc(c)*&
                   (uw%n(c)*phi_x(c)+vw%n(c)*phi_y(c)+(ww%n(c))*phi_z(c)) 
      end do
    end if
    call GraPhi(VAR1x,1,VAR2x,.TRUE.)
    call GraPhi(VAR1y,2,VAR2y,.TRUE.)
    call GraPhi(VAR1z,3,VAR2z,.TRUE.)

    do c=1,NC
      b(c) = b(c) + (VAR2x(c)+VAR2y(c)+VAR2z(c))*volume(c)
    end do

    !------------------------------------------------------------------!
    !   Here we clean up transport equation from the false diffusion   !
    !------------------------------------------------------------------!
    if(SIMULA==EBM.and.MODE/=HYB) then
      do s=1,NS

        c1=SideC(1,s)
        c2=SideC(2,s)

        VISeff = (fF(s)*VISt(c1)+(1.0-fF(s))*VISt(c2)) 

        phi_xS = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2)
        phi_yS = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
        phi_zS = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)
        Fex=VISeff*(phi_xS*Sx(s) + phi_yS*Sy(s) + phi_zS*Sz(s))
        A0 = VISeff * Scoef(s)
        Fim = (   phi_xS*Dx(s)      &
                + phi_yS*Dy(s)      &
                + phi_zS*Dz(s))*A0

        b(c1) = b(c1) - VISeff*(phi%n(c2)-phi%n(c1))*Scoef(s)- Fex + Fim
        if(c2  > 0) then
          b(c2) = b(c2) + VISeff*(phi%n(c2)-phi%n(c1))*Scoef(s)+ Fex - Fim
        end if
      end do
    end if
  end if

  !---------------------------------!
  !     Temporal discretization     !
  !---------------------------------!

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
  !  is handled via the linear system of equations 

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

  ! Two time levels; Linear interpolation
  if(INERT == LIN) then
    do c=1,NC
      A0 = DENc(material(c))*volume(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(INERT == PAR) then
    do c=1,NC
      A0 = DENc(material(c))*volume(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * A0
      b(c) = b(c) + 2.0*A0 * phi % o(c) - 0.5*A0 * phi % oo(c)
    end do
  end if

  !  do c=1,NC                                         ! 2mat
  !    if(StateMat(material(c))==SOLID) then           ! 2mat
  !      b(c)          = 0.0                           ! 2mat
  !    end if                                          ! 2mat
  !  end do                                            ! 2mat

  !---- Disconnect the SOLID cells from FLUID system   ! 2mat
  !  do s=1,NS                                         ! 2mat
  !    c1=SideC(1,s)                                   ! 2mat
  !    c2=SideC(2,s)                                   ! 2mat
  !    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat
  !      if(c2 > 0) then ! => not buffer               ! 2mat
  !        if(StateMat(material(c1)) == SOLID) then    ! 2mat
  !          A % val(A % pos(1,s)) = 0.0               ! 2mat 
  !          A % val(A % pos(2,s)) = 0.0               ! 2mat
  !        end if                                      ! 2mat 
  !        if(StateMat(material(c2)) == SOLID) then    ! 2mat
  !          A % val(A % pos(2,s)) = 0.0               ! 2mat
  !          A % val(A % pos(1,s)) = 0.0               ! 2mat 
  !        end if                                      ! 2mat 
  !      else            ! => buffer region            ! 2mat 
  !        if(StateMat(material(c1)) == SOLID .or.  &  ! 2mat
  !           StateMat(material(c2)) == SOLID) then    ! 2mat
  !          A % bou(c2) = 0.0                         ! 2mat
  !        end if                                      ! 2mat 
  !      end if                                        ! 2mat
  !    end if                                          ! 2mat
  !  end do                                            ! 2mat

   if(SIMULA==EBM) then 
     call SourcesEBM(var)
   else
     call SourcesHJ(var)        
   end if                

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

  if(ALGOR == SIMPLE)   miter=30
  if(ALGOR == FRACT)    miter=5

  niter=miter
  call cg(NC, Nbc, A,           &
           phi % n, b, PREC,    &
           niter,phi % STol, res(var), error)

  !k-eps  if(var == 1) then
  !k-eps    write(LineRes(17:28), '(1PE12.3)') res(var) 
  !k-eps    write(LineRes(77:80), '(I4)')      niter 
  !k-eps  end if

  if(this_proc < 2) write(*,*) 'Var ', var, res(var), niter 

  if(var == 13) then
    do c= 1, NC
      phi%n(c) = phi%n(c) 
     if( phi%n(c)<0.0)then
       phi%n(c) = phi%o(c)
     end if
    end do
  end if

  if(var == 6.or.var == 7.or.var == 8) then
    do c= 1, NC
      phi % n(c) = phi % n(c) 
      if(phi % n(c) < 0.0) then
        phi % n(c) = phi % o(c)
      end if
    end do
  end if

  call Exchng(phi % n)

  end subroutine
