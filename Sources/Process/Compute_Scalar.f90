!==============================================================================!
  subroutine Compute_Scalar(grid, var, phi)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod
  use Grid_Mod
  use Var_Mod
  use Work_Mod, only: phi_x => r_cell_01,  &
                      phi_y => r_cell_02,  &
                      phi_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Grid_Type) :: grid
  integer         :: var
  type(Var_Type)  :: phi
!----------------------------------[Calling]-----------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------! 
  integer :: n, c, s, c1, c2, niter, miter, mat
  real    :: A0, A12, A21, error
  real    :: CONeff1, FUex1, FUim1, phixS1, phiyS1, phizS1
  real    :: CONeff2, FUex2, FUim2, phixS2, phiyS2, phizS2
  real    :: Stot, phis, CAPs, Prt1, Prt2
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

!  if(SIMULA==EBM) then
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

  ! Old values (o and oo)
  if(ini.lt.2) then
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
    if(BLEND_TEM(mat) /= NO) then
      call Compute_Minimum_Maximum(grid, phi % n)  ! or phi % o ???
      goto 1
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

 
    if(BLEND_TEM(material(c1)) /= NO .or.  &
       BLEND_TEM(material(c2)) /= NO) then
      call Advection_Scheme(grid, phis, s, phi % n,     &
                      phi_x, phi_y, phi_z,              &
                      grid % dx, grid % dy, grid % dz,  &
                      max(BLEND_TEM(material(c1)),      &
                          BLEND_TEM(material(c2))) ) 
    end if 

    CAPs = grid % f(s)       * CAPc(material(c1)) &
         + (1.0-grid % f(s)) * CAPc(material(c2))

    ! Central differencing for advection
    if(ini.eq.1) then
      if(c2.gt.0) then
        phi % a_o(c1) = phi % a_o(c1)-Flux(s)*phis*CAPs
        phi % a_o(c2) = phi % a_o(c2)+Flux(s)*phis*CAPs
      else
        phi % a_o(c1)=phi % a_o(c1)-Flux(s)*phis*CAPs
      endif
    end if
    if(c2.gt.0) then
      phi % a(c1) = phi % a(c1)-Flux(s)*phis*CAPs
      phi % a(c2) = phi % a(c2)+Flux(s)*phis*CAPs
    else
      phi % a(c1) = phi % a(c1)-Flux(s)*phis*CAPs
    endif
 
    ! Upwind
    if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
      if(Flux(s).lt.0) then   ! from c2 to c1
        phi % c(c1) = phi % c(c1)-Flux(s)*phi % n(c2) * CAPs
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+Flux(s)*phi % n(c2) * CAPs
        endif
      else
        phi % c(c1) = phi % c(c1)-Flux(s)*phi % n(c1) * CAPs
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+Flux(s)*phi % n(c1) * CAPs
        endif
      end if
    end if   ! BLEND_TEM
  end do  ! through sides

  !-----------------------------!
  !                             !
  !   Temporal discretization   !
  !                             !
  !-----------------------------!
   
  ! Adams-Bashforth scheeme for convective fluxes
  if(CONVEC.eq.AB) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tem(material(c)) * &
                    (1.5*phi % a_o(c) - 0.5*phi % a_oo(c) - phi % c(c))
    end do
  endif

  ! Crank-Nicholson scheeme for convective fluxes
  if(CONVEC.eq.CN) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tem(material(c)) * &
                    (0.5 * ( phi % a(c) + phi % a_o(c) ) - phi % c(c))
    end do
  endif

  ! Fully implicit treatment of convective fluxes
  if(CONVEC.eq.FI) then
    do c = 1, grid % n_cells
      b(c) = b(c) + URFC_Tem(material(c)) * &
                    (phi % a(c)-phi % c(c))
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
     
    if(SIMULA/=LES.or.SIMULA/=DNS) then
      Prt1 = 1.0/( 0.5882 + 0.228*(VISt(c1)/(VISc+1.0e-12)) - 0.0441*                  &
            (VISt(c1)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*( VISc/(VISt(c1)+1.0e-12) ))) )
      Prt2 = 1.0/( 0.5882 + 0.228*(VISt(c2)/(VISc+1.0e-12)) - 0.0441*                  &
           (VISt(c2)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*( VISc/(VISt(c2)+1.0e-12) ))) )
      Prt = fF(s)*Prt1 + (1.0-fF(s))*Prt2
    end if

    ! Gradients on the cell face 
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      if(material(c1) == material(c2)) then
        phixS1 = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2) 
        phiyS1 = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
        phizS1 = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)
        phixS2 = phixS1 
        phiyS2 = phiyS1 
        phizS2 = phizS1 
        CONeff1 =      grid % f(s) * ( CONc(material(c1))                 &
                                     + CAPc(material(c1))*VISt(c1)/Prt )  &
                + (1.-grid % f(s)) * ( CONc(material(c2))                 &
                                     + CAPc(material(c2))*VISt(c2)/Prt )
        CONeff2 = CONeff1 
      else 
        phixS1 = phi_x(c1) 
        phiyS1 = phi_y(c1) 
        phizS1 = phi_z(c1) 
        phixS2 = phi_x(c2) 
        phiyS2 = phi_y(c2) 
        phizS2 = phi_z(c2) 
        CONeff1 =   CONc(material(c1))                 &
                  + CAPc(material(c1))*VISt(c1)/Prt   
        CONeff2 =   CONc(material(c2))                 &
                  + CAPc(material(c2))*VISt(c2)/Prt   
      end if
    else
      phixS1 = phi_x(c1) 
      phiyS1 = phi_y(c1) 
      phizS1 = phi_z(c1) 
      phixS2 = phixS1 
      phiyS2 = phiyS1 
      phizS2 = phizS1 
      CONeff1 =   CONc(material(c1))                 &
                + CAPc(material(c1))*VISt(c1)/Prt   
      CONeff2 = CONeff1 
    endif


    if(SIMULA == ZETA      .or.  &
       SIMULA == K_EPS_VV  .or.  &
       SIMULA == K_EPS     .or.  &
       SIMULA == HYB_ZETA) then
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
          CONeff1 = CONwall(c1)
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
    FUim1 = CONeff1*Scoef(s)*    &
            (  phixS1*grid % dx(s)      &
             + phiyS1*grid % dy(s)      &
             + phizS1*grid % dz(s) )
    FUim2 = CONeff2*Scoef(s)*    &
            (  phixS2*grid % dx(s)      &
             + phiyS2*grid % dy(s)      &
             + phizS2*grid % dz(s) )

    ! Straight diffusion part 
    if(ini.lt.2) then
      if(c2.gt.0) then
        if(material(c1) == material(c2)) then
          phi % d_o(c1) = phi % d_o(c1)  &
                       + CONeff1*Scoef(s)*(phi % n(c2) - phi % n(c1)) 
          phi % d_o(c2) = phi % d_o(c2)  &
                       - CONeff2*Scoef(s)*(phi % n(c2) - phi % n(c1))   
        else
          phi % d_o(c1) = phi % d_o(c1)            &
                       + 2.*CONc(material(c1))   &
                       * Scoef(s) * (phiside(s) - phi % n(c1)) 
          phi % d_o(c2) = phi % d_o(c2)            &
                       - 2.*CONc(material(c2))   &
                       * Scoef(s)*(phi % n(c2) - phiside(s))   
        end if
      else
        if(TypeBC(c2).ne.SYMMETRY) then 
          if(material(c1) == material(c2)) then
            phi % d_o(c1) = phi % d_o(c1)  &
                + CONeff1*Scoef(s)*(phi % n(c2) - phi % n(c1))   
          else
            phi % d_o(c1) = phi % d_o(c1)  &
                + 2.*CONc(material(c1))*Scoef(s)*(phiside(s)-phi % n(c1)) 
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
    if( (DIFFUS.eq.CN) .or. (DIFFUS.eq.FI) ) then

      if(DIFFUS .eq. CN) then       ! Crank Nicholson
        if(material(c1) == material(c2)) then
          A12 = .5*CONeff1*Scoef(s)  
          A21 = .5*CONeff2*Scoef(s)  
        else
          A12 = CONc(material(c1))*Scoef(s)  
          A21 = CONc(material(c2))*Scoef(s)  
        end if
      end if

      if(DIFFUS .eq. FI) then       ! Fully implicit
        if(material(c1) == material(c2)) then
          A12 = CONeff1*Scoef(s)  
          A21 = CONeff2*Scoef(s)  
        else
          A12 = 2.*CONc(material(c1))*Scoef(s)  
          A21 = 2.*CONc(material(c2))*Scoef(s)  
        end if
      end if

      if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
        A12 = A12  - min(Flux(s), 0.0)
        A21 = A21  + max(Flux(s), 0.0)
      endif
                
      ! Fill the system matrix
      if(c2.gt.0) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + A21
        if(material(c1) == material(c2)) then
          A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A12
          A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A21
        else
          b(c1) = b(c1) + A12*phiside(s)
          b(c2) = b(c2) + A21*phiside(s)
        end if
      else if(c2.lt.0) then

        ! Outflow is included because of the flux 
        ! corrections which also affects velocities
        if( (TypeBC(c2).eq.INFLOW).or.    &
            (TypeBC(c2).eq.WALL).or.      &
            (TypeBC(c2).eq.CONVECT) ) then    
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          b(c1)  = b(c1)  + A12 * phi % n(c2)

        ! Buffer: System matrix and parts belonging 
        ! to other subdomains are filled here.
        else if(TypeBC(c2).eq.BUFFER) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + A12
          if(material(c1) == material(c2)) then
            A % bou(c2) = -A12  ! cool parallel stuff
          else
            b(c1) = b(c1) + A12*phiside(s)
          end if

        ! In case of wallflux 
        else if(TypeBC(c2).eq.WALLFL) then
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
  if(DIFFUS.eq.AB) then 
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % d_o(c) - 0.5*phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(DIFFUS.eq.CN) then 
    do c = 1, grid % n_cells
      b(c)  = b(c) + 0.5*phi % d_o(c)
    end do  
  end if

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(CROSS.eq.AB) then
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % c_o(c) - 0.5*phi % c_oo(c)
    end do 
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(CROSS.eq.CN) then
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
  if(CROSS.eq.FI) then
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
  if(INERT.eq.LIN) then
    do c = 1, grid % n_cells
      A0 = CAPc(material(c)) * DENc(material(c)) * grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + A0
      b(c)  = b(c) + A0*phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(INERT.eq.PAR) then
    do c = 1, grid % n_cells
      A0 = CAPc(material(c)) * DENc(material(c)) * grid % vol(c)/dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * A0
      b(c)  = b(c) + 2.0 * A0 * phi % o(c) - 0.5 * A0 * phi % oo(c)
    end do
  end if

  if(SIMULA==EBM.or.SIMULA==HJ) then
    if(MODE/=HYB) then
      do c = 1, grid % n_cells
        VAR1x(c) = -0.22*Tsc(c) *&
                   (uu%n(c)*phi_x(c)+uv%n(c)*phi_y(c)+uw%n(c)*phi_z(c))
        VAR1y(c) = -0.22*Tsc(c)*&
                   (uv%n(c)*phi_x(c)+vv%n(c)*phi_y(c)+vw%n(c)*phi_z(c))
        VAR1z(c) = -0.22*Tsc(c)*&
                   (uw%n(c)*phi_x(c)+vw%n(c)*phi_y(c)+ww%n(c)*phi_z(c))
      end do
      call GraPhi(grid, VAR1x,1,VAR2x,.TRUE.)
      call GraPhi(grid, VAR1y,2,VAR2y,.TRUE.)
      call GraPhi(grid, VAR1z,3,VAR2z,.TRUE.)
      do c = 1, grid % n_cells
        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
      end do

      !------------------------------------------------------------------!
      !   Here we clean up transport equation from the false diffusion   !
      !------------------------------------------------------------------!
      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        Prt1 = 1.0/( 0.5882 + 0.228*(VISt(c1)/(VISc+1.0e-12)) - 0.0441*                  &
              (VISt(c1)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*( VISc/(VISt(c1)+1.0e-12) ))) )
        Prt2 = 1.0/( 0.5882 + 0.228*(VISt(c2)/(VISc+1.0e-12)) - 0.0441*                  &
              (VISt(c2)/(VISc+1.0e-12))**2.0*(1.0 - exp(-5.165*( VISc/(VISt(c2)+1.0e-12) ))) )

        Prt = fF(s)*Prt1 + (1.0-fF(s))*Prt2
        if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
          phixS1 = fF(s)*phi_x(c1) + (1.0-fF(s))*phi_x(c2) 
          phiyS1 = fF(s)*phi_y(c1) + (1.0-fF(s))*phi_y(c2)
          phizS1 = fF(s)*phi_z(c1) + (1.0-fF(s))*phi_z(c2)
          phixS2 = phixS1 
          phiyS2 = phiyS1 
          phizS2 = phizS1 
          CONeff1 =       grid % f(s)  * (CAPc(material(c1))*VISt(c1)/Prt )  &
                  + (1. - grid % f(s)) * (CAPc(material(c2))*VISt(c2)/Prt )
          CONeff2 = CONeff1 
        else
          phixS1 = phi_x(c1) 
          phiyS1 = phi_y(c1) 
          phizS1 = phi_z(c1) 
          phixS2 = phixS1 
          phiyS2 = phiyS1 
          phizS2 = phizS1 
          CONeff1 = CAPc(material(c1))*VISt(c1)/Prt   
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
        FUim1 = CONeff1*Scoef(s)*    &
                (  phixS1*grid % dx(s)      &
                 + phiyS1*grid % dy(s)      &
                 + phizS1*grid % dz(s) )
        FUim2 = CONeff2*Scoef(s)*    &
                (  phixS2*grid % dx(s)      &
                 + phiyS2*grid % dy(s)      &
                 + phizS2*grid % dz(s) )

        b(c1) = b(c1) - CONeff1*(phi%n(c2)-phi%n(c1))*Scoef(s)- FUex1 + FUim1
        if(c2  > 0) then
         b(c2) = b(c2) + CONeff1*(phi%n(c2)-phi%n(c1))*Scoef(s)+ FUex2 - FUim2
        end if
      end do
    end if  
  end if  

  call User_Source(grid)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - phi % URF) * phi % n(c)  &
                                      / phi % URF
    A % val(A % dia(c)) = A % val(A % dia(c)) / phi % URF
  end do  

  if(ALGOR == SIMPLE)   miter=10
  if(ALGOR == FRACT)    miter=5 

  niter=miter
  call bicg(A, phi % n, b,            &
            PREC, niter, phi % STol,  &
            res(var), error)
  write(LineRes(65:76),  '(1PE12.3)') res(var)
  write(LineRes(93:96),  '(I4)')      niter       

  call Exchange(grid, phi % n)

  end subroutine
