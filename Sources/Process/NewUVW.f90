!==============================================================================!
  subroutine NewUVW(grid, var, ui,                                  &
                    ui_i, ui_j, ui_k,                            &
                    Si, Sj, Sk,                                     &
                    Di, Dj, Dk,                                     &
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
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
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
  integer :: s, c, c1, c2, niter, miter, mat
  real    :: Fex, Fim 
  real    :: uis
  real    :: A0, A12, A21
  real    :: error
  real    :: VISeff, VIStS, Fstress 
  real    :: ui_iS,ui_jS,ui_kS,uj_iS,uk_iS
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

  ! Old values (o) and older than old (oo)
  if(ini == 1) then
    do c=1,grid % n_cells
      ui % oo(c)  = ui % o(c)
      ui % o (c)  = ui % n(c)
      ui % Coo(c) = ui % Co(c)
      ui % Co (c) = 0.0 
      ui % Doo(c) = ui % Do(c)
      ui % Do (c) = 0.0 
      ui % Xoo(c) = ui % Xo(c)
      ui % Xo (c) = ui % X(c) 
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Compute phimax and phimin
  do mat=1,grid % n_materials
    if(BLEND(mat) /= NO) then
      call Compute_Minimum_Maximum(grid, ui % n)  ! or ui % o ???
      goto 1
    end if
  end do

  ! New values
1 do c=1,grid % n_cells
    ui % C(c)    = 0.0
    ui % X(c)    = 0.0
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s=1,grid % n_faces

    c1=SideC(1,s)
    c2=SideC(2,s) 

    ! Central differencing
    uis=f(s)*ui % n(c1) + (1.0-f(s))*ui % n(c2)

    if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
      call Advection_Scheme(grid, uis, s, ui % n, ui_i, ui_j, ui_k, Di, Dj, Dk, &
                           max(BLEND(material(c1)),BLEND(material(c2))) ) 
    end if 
    
    ! Central differencing for advection
    if(ini == 1) then 
      if(c2  > 0) then
        ui % Co(c1)=ui % Co(c1)-Flux(s)*uis
        ui % Co(c2)=ui % Co(c2)+Flux(s)*uis
      else
        ui % Co(c1)=ui % Co(c1)-Flux(s)*uis
      endif 
    end if

    if(c2  > 0) then
      ui % C(c1)=ui % C(c1)-Flux(s)*uis
      ui % C(c2)=ui % C(c2)+Flux(s)*uis
    else
      ui % C(c1)=ui % C(c1)-Flux(s)*uis
    endif 

    ! Upwind 
    if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
      if(Flux(s)  < 0) then   ! from c2 to c1
        ui % X(c1)=ui % X(c1)-Flux(s)*ui % n(c2)
        if(c2  > 0) then
          ui % X(c2)=ui % X(c2)+Flux(s)*ui % n(c2)
        endif
      else 
        ui % X(c1)=ui % X(c1)-Flux(s)*ui % n(c1)
        if(c2  > 0) then
          ui % X(c2)=ui % X(c2)+Flux(s)*ui % n(c1)
        endif
      end if
    end if   ! BLEND 
  end do    ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(CONVEC == AB) then
    do c=1,grid % n_cells
      b(c) = b(c) + URFC(material(c)) * & 
                    (1.5*ui % Co(c) - 0.5*ui % Coo(c) - ui % X(c))
    end do  
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(CONVEC == CN) then
    do c=1,grid % n_cells
      b(c) = b(c) + URFC(material(c)) * & 
                    (0.5 * ( ui % C(c) + ui % Co(c) ) - ui % X(c))
    end do  
  endif

  ! Fully implicit treatment of convective fluxes 
  if(CONVEC == FI) then
    do c=1,grid % n_cells
      b(c) = b(c) + URFC(material(c)) * & 
                    (ui % C(c) - ui % X(c))
    end do  
  end if     
          
  ! New values
  do c=1,grid % n_cells
    ui % X(c) = 0.0
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

    c1=SideC(1,s)
    c2=SideC(2,s)   

    VISeff = fF(s)*VISt(c1)+(1.0-fF(s))*VISt(c2) + VISc

    if(SIMULA==HYB_ZETA) VISeff = fF(s)*VISt_eff(c1)+(1.0-fF(s))*VISt_eff(c2) + VISc

    if(c2 < 0 .and. SIMULA == LES) then
      if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
        VISeff = VISwall(c1)
      end if
    end if 

    if(SIMULA == ZETA.or.SIMULA==K_EPS_VV.or.(SIMULA == K_EPS.and.MODE==HRe).or.SIMULA==HYB_ZETA) then
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
          VISeff = VISwall(c1)
        end if
      end if
    end if

    ! Add influence of Re stresses for EBM
    if(SIMULA == EBM.or.SIMULA == HJ) then
      if(MODE /= HYB) then        
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
        ui % Do(c1) = ui % Do(c1) + (ui % n(c2)-ui % n(c1))*A0   
        ui % Do(c2) = ui % Do(c2) - (ui % n(c2)-ui % n(c1))*A0    
      else
        if(TypeBC(c2) /= SYMMETRY) then
          ui % Do(c1) = ui % Do(c1) + (ui % n(c2)-ui % n(c1))*A0   
        end if 
      end if 
    end if

    ! Cross diffusion part
    ui % X(c1) = ui % X(c1) + Fex - Fim + Fstress
    if(c2  > 0) then
      ui % X(c2) = ui % X(c2) - Fex + Fim - Fstress
    end if 

     ! Compute the coefficients for the sysytem matrix
    if( (DIFFUS == CN) .or. (DIFFUS == FI) ) then  
      if(DIFFUS  ==  CN) then       ! Crank Nicholson
        A12 = 0.5 * A0 
        A21 = 0.5 * A0 
      end if

      if(DIFFUS  ==  FI) then       ! Fully implicit
        A12 = A0 
        A21 = A0
      end if
    
      if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
        A12 = A12  - min(Flux(s), real(0.0)) 
        A21 = A21  + max(Flux(s), real(0.0))
      endif

      ! Fill the system matrix
      if(c2  > 0) then
        A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A12
        A % val(A % dia(c1))    = A % val(A % dia(c1))    + A12
        A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A21
        A % val(A % dia(c2))    = A % val(A % dia(c2))    + A21
      else if(c2  < 0) then

        ! Outflow is not included because it was causing problems     
        if((TypeBC(c2) == INFLOW).or.                                &
           (TypeBC(c2) == WALL).or.                                  &
           (TypeBC(c2) == CONVECT).or.                               &
           (TypeBC(c2) == WALLFL)) then                                
           ! (TypeBC(c2) == OUTFLOW) ) then   
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1) = b(c1) + A12 * ui % n(c2)
        else if(TypeBC(c2) == BUFFER) then  
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          A % bou(c2) = -A12  ! cool parallel stuff
        endif
      end if     
    end if
  end do  ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  !---------------------------------!
  !
  ! Add Re stress influence on momentum
  ! This is an alternative way to implement RSM. 
  !
  !  if(SIMULA == EBM.or.SIMULA == HJ) then
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
  if(SIMULA == EBM.or.SIMULA == HJ) then
    if(MODE /= HYB) then
      do s=1,grid % n_faces
        c1=SideC(1,s)
        c2=SideC(2,s)

        VIStS = (fF(s)*VISt(c1)+(1.0-fF(s))*VISt(c2))
        A0 = Scoef(s)*VIStS 
        VISeff = VIStS

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
  if(DIFFUS == AB) then 
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * ui % Do(c) - 0.5 * ui % Doo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * ui % Do(c)
    end do  
  end if
             
  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c=1,grid % n_cells
      b(c) = b(c) + 1.5 * ui % Xo(c) - 0.5 * ui % Xoo(c)
    end do 
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c=1,grid % n_cells
      b(c) = b(c) + 0.5 * ui % X(c) + 0.5 * ui % Xo(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
    do c=1,grid % n_cells
      b(c) = b(c) + ui % X(c)
    end do 
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; linear interpolation
  if(INERT == LIN) then
    do c=1,grid % n_cells
      A0 = DENc(material(c))*grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c) = b(c) + A0 * ui % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(INERT == PAR) then
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
      b(c) = b(c)  + PdropX(material(c)) * grid % vol(c)
    end do
  else if(ui % name == 'V') then
    do c=1,grid % n_cells
      b(c) = b(c)  + PdropY(material(c)) * grid % vol(c)
    end do
  else if(ui % name == 'W') then
    do c=1,grid % n_cells
      b(c) = b(c)  + PdropZ(material(c)) * grid % vol(c)
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
  if(HOT == YES) call User_Force(grid, ui, a, b)

!  do c=1,grid % n_cells                                         ! 2mat
!    if(StateMat(material(c))==SOLID) then           ! 2mat
!      b(c)          = 0.0                           ! 2mat
!    end if                                          ! 2mat
!  end do                                            ! 2mat

!---- Disconnect the SOLID cells from FLUID system  ! 2mat
!  do s=1,grid % n_faces                             ! 2mat
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

  !-----------------------------------!
  !                                   !
  !   Solve the equations for u,v,w   !
  !                                   !    
  !-----------------------------------!
  do c=1,grid % n_cells
    A % sav(c) = A % val(A % dia(c))
    b(c) = b(c) + A % val(A % dia(c)) * (1.0-U % URF)*ui % n(c) / U % URF
    A % val(A % dia(c)) = A % val(A % dia(c)) / U % URF
  end do

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5

  niter=miter

  call cg(grid % n_cells, grid % n_bnd_cells, A,           & 
          ui % n, b, PREC,        &
          niter,U % STol, res(var), error)

  if(ui % name == 'U') then
    write(LineRes(17:28), '(1PE12.3)') res(var) 
    write(LineRes(77:80), '(I4)')      niter 
  end if
  if(ui % name == 'V') then
    write(LineRes(29:40), '(1PE12.3)') res(var) 
    write(LineRes(81:84), '(I4)')      niter 
  end if
  if(ui % name == 'W') then
    write(LineRes(41:52), '(1PE12.3)') res(var) 
    write(LineRes(85:88), '(I4)')      niter 
  end if

  call Exchng(ui % n)

  end subroutine
