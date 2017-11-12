!==============================================================================!
  subroutine Compute_Pressure_Simple(grid)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the S.I.M.P.L.E. method.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, niter
  real    :: Us, Vs, Ws, DENs, A12      
  real    :: Pmax, Pmin
  real    :: error
  real    :: SMDPN
  real    :: dPxi, dPyi, dPzi
!------------------------------------------------------------------------------!
!     
!  The form of equations which I am solving:    
!     
!     /               /            
!    |               |             
!    | rho u dS = dt | GRAD pp dS
!    |               |             
!   /               /              
!
!  Dimension of the system under consideration
!   
!     [App] {pp} = {bpp}               [kg/s]
!   
!  Dimensions of certain variables
!
!     APP            [ms]
!     PP             [kg/ms^2]
!     b              [kg/s]
!     Flux           [kg/s]
!   
!==============================================================================!

  ! Initialize matrix and right hand side
  A % val = 0.0
  b = 0.0 

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !-----------------------------------------!
  PP % n = 0.0 

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s=1, NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    DENs =      f(s)  * DENc(material(c1))     &
         + (1.0-f(s)) * DENc(material(c2))

    ! Handle materials
    ! if( StateMat(material(c1))==FLUID .and.        &
    !   StateMat(material(c2))==FLUID) then
    !   DENs =      f(s)  * DENc(material(c1))       &
    !        + (1.0-f(s)) * DENc(material(c2))
    !   else if( StateMat(material(c1))==FLUID .and. &
    !            StateMat(material(c2))==SOLID) then
    !     DENs = DENc(material(c1)) 
    !   else if( StateMat(material(c1))==SOLID .and. &
    !            StateMat(material(c2))==FLUID) then
    !     DENs = DENc(material(c2)) 
    !   else
    !     DENs =      f(s)  * DENc(material(c1))     &
    !          + (1.0-f(s)) * DENc(material(c2))
    ! end if  

    ! Face is inside the domain
    if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then

      SMDPN = (  grid % sx(s)*grid % sx(s)   &
               + grid % sy(s)*grid % sy(s)   &
               + grid % sz(s)*grid % sz(s) ) &
            / (  grid % sx(s)*grid % dx(s)   & 
               + grid % sy(s)*grid % dy(s)   &
               + grid % sz(s)*grid % dz(s) )  

      ! Interpolate velocity 
      Us = f(s) * U % n(c1) + (1.0-f(s)) * U % n(c2)
      Vs = f(s) * V % n(c1) + (1.0-f(s)) * V % n(c2)
      Ws = f(s) * W % n(c1) + (1.0-f(s)) * W % n(c2)

      ! Calculate coeficients for the system matrix
      if(c2  > 0) then 
        A12 = 0.5 * DENs * SMDPN *           &
           (  grid % vol(c1) / A % sav(c1)       &
            + grid % vol(c2) / A % sav(c2) )  
        A % val(A % pos(1,s))  = -A12
        A % val(A % pos(2,s))  = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
        A % val(A % dia(c2)) = A % val(A % dia(c2)) +  A12
      else 
        A12 = 0.5 * DENs * SMDPN *           &
             (  grid % vol(c1) / A % sav(c1)     &
              + grid % vol(c2) / A % sav(c2) )  
        A % bou(c2)  = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      end if 

      ! Interpolate pressure gradients
      dPxi=.5*( Px(c1) + Px(c2) )*grid % dx(s)
      dPyi=.5*( Py(c1) + Py(c2) )*grid % dy(s)
      dPzi=.5*( Pz(c1) + Pz(c2) )*grid % dz(s)

      ! Calculate flux through cell face
      Flux(s) = DENs * (  Us*grid % sx(s)       &
                        + Vs*grid % sy(s)       &
                        + Ws*grid % sz(s) )     &
              + A12 * (P % n(c1) - P % n(c2))   &
              + A12 * (dPxi + dPyi + dPzi)                            

      b(c1)=b(c1)-Flux(s)
      if(c2  > 0) b(c2)=b(c2)+Flux(s)

    ! Side is on the boundary
    else ! (c2 < 0)

      if(TypeBC(c2) == INFLOW) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        Flux(s) = DENs * (  Us * grid % sx(s)  &
                          + Vs * grid % sy(s)  &
                          + Ws * grid % sz(s) )
        b(c1) = b(c1)-Flux(s)
      else if(TypeBC(c2) == OUTFLOW .or.   &
              TypeBC(c2) == CONVECT) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        Flux(s) = DENs * (  Us*grid % sx(s)  &
                          + Vs*grid % sy(s)  &
                          + Ws*grid % sz(s) )
        b(c1) = b(c1)-Flux(s)
        SMDPN = (  grid % sx(s) * grid % sx(s)   &
                 + grid % sy(s) * grid % sy(s)   &
                 + grid % sz(s) * grid % sz(s) ) &
              / (  grid % sx(s) * grid % dx(s)   &
                 + grid % sy(s) * grid % dy(s)   &
                 + grid % sz(s) * grid % dz(s) )  
        A12 = DENs * SMDPN * grid % vol(c1) / A % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      else if(TypeBC(c2) == PRESSURE) then
        Us = U % n(c1)
        Vs = V % n(c1)
        Ws = W % n(c1)
        Flux(s) = DENs * (  Us * grid % sx(s)  &
                          + Vs * grid % sy(s)  &
                          + Ws * grid % sz(s) )
        b(c1) = b(c1)-Flux(s)
        SMDPN = ( grid % sx(s) * grid % sx(s)   &
                + grid % sy(s) * grid % sy(s)   &
                + grid % sz(s) * grid % sz(s) ) &
              / ( grid % sx(s) * grid % dx(s)   &
                + grid % sy(s) * grid % dy(s)   &
                + grid % sz(s) * grid % dz(s) )
        A12 = DENs * SMDPN * grid % vol(c1) / A % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      else  ! it is SYMMETRY
        Flux(s) = 0.0
      end if
    end if

1 end do

  errmax=0.0
  do c=1,NC
    errmax=max(errmax, abs(b(c)))
  end do

  ! Don't solve the pressure corection too accurate.
  ! Value 1.e-18 blows the solution.
  ! Value 1.e-12 keeps the solution stable
  niter=40
  call bicg(NC, NbC, A,                     &
        PP % n, b, PREC, niter, PP % STol,  &
        res(4), error) 
  write(LineRes(53:64),  '(1PE12.3)') res(4)
  write(LineRes(89:92),  '(I4)')      niter

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = 1, NC
    P % n(c)  =  P % n(c)  +  P % URF  *  PP % n(c)
  end do

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  !Pmax  = maxval(P % n(1:NC))
  !Pmin  = minval(P % n(1:NC))
 
  !call glomax(Pmax)
  !call glomin(Pmin)
 
  !P % n = P % n - 0.5*(Pmax+Pmin)
 
  call Exchng(PP % n)

  end subroutine
