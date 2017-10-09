!======================================================================!
  subroutine CalcPS
!----------------------------------------------------------------------!
!   Forms and solves pressure equation for the S.I.M.P.L.E. method.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: s, c, c1, c2, niter
  real    :: Us, Vs, Ws, DENs, A12      
  real    :: Pmax, Pmin
  real    :: error
  real    :: SMDPN
  real    :: dPxi, dPyi, dPzi
!----------------------------------------------------------------------!
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
!======================================================================!

  A % val = 0.0

!------------------------------------------!
!      Initialize the source term for      !
!     the pressure correction equation     !
!------------------------------------------!
  b = 0.0 

!---------------------------------------------!
!     Initialize the pressure corrections     !
!---------------------------------------------!
  PP % n = 0.0 

!-----------------------------------------------------!
!     Calculate the mass fluxes on the cell faces     !
!-----------------------------------------------------!
  do s=1, NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    DENs =      f(s)  * DENc(material(c1))     &
         + (1.0-f(s)) * DENc(material(c2))
!---- handle materials
!    if( StateMat(material(c1))==FLUID .and.      &
!        StateMat(material(c2))==FLUID) then
!      DENs =      f(s)  * DENc(material(c1))     &
!           + (1.0-f(s)) * DENc(material(c2))
!    else if( StateMat(material(c1))==FLUID .and. &
!             StateMat(material(c2))==SOLID) then
!      DENs = DENc(material(c1)) 
!    else if( StateMat(material(c1))==SOLID .and. &
!             StateMat(material(c2))==FLUID) then
!      DENs = DENc(material(c2)) 
!    else
!      DENs =      f(s)  * DENc(material(c1))     &
!           + (1.0-f(s)) * DENc(material(c2))
!    end if  

!---- side is inside the domain
    if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then

      SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
            / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )  

!---- interpoliraj gustocu i brzine
      Us = f(s) * U % n(c1) + (1.0-f(s)) * U % n(c2)
      Vs = f(s) * V % n(c1) + (1.0-f(s)) * V % n(c2)
      Ws = f(s) * W % n(c1) + (1.0-f(s)) * W % n(c2)

!---- calculate coeficients for the system matrix
      if(c2  > 0) then 
        A12 = 0.5 * DENs * SMDPN *                                  &
           (  volume(c1) / A % sav(c1)                              &
            + volume(c2) / A % sav(c2) )  
        A % val(A % pos(1,s))  = -A12
        A % val(A % pos(2,s))  = -A12
       A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
        A % val(A % dia(c2)) = A % val(A % dia(c2)) +  A12
      else 
        A12 = 0.5 * DENs * SMDPN *                                  &
             (  volume(c1) / A % sav(c1)                              &
              + volume(c2) / A % sav(c2) )  
        A % bou(c2)  = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      end if 

!---- interpoliraj razliku pritiska
      dPxi=.5*( Px(c1) + Px(c2) )*Dx(s)
      dPyi=.5*( Py(c1) + Py(c2) )*Dy(s)
      dPzi=.5*( Pz(c1) + Pz(c2) )*Dz(s)

!---- now calculate the flux through cell face
      Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )           &
            + A12 * (P % n(c1) - P % n(c2))                       &
            + A12 * (dPxi+dPyi+dPzi)                            

      b(c1)=b(c1)-Flux(s)
      if(c2  > 0) b(c2)=b(c2)+Flux(s)

!----- side is on the boundary
    else ! (c2 < 0)

      if(TypeBC(c2) == INFLOW) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )
        b(c1) = b(c1)-Flux(s)
      else if(TypeBC(c2) == OUTFLOW .or.   &
              TypeBC(c2) == CONVECT) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )
        b(c1) = b(c1)-Flux(s)
        SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
              / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )  
        A12 = DENs * SMDPN * volume(c1) / A % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      else if(TypeBC(c2) == PRESSURE) then
        Us = U % n(c1)
        Vs = V % n(c1)
        Ws = W % n(c1)
        Flux(s) = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )
        b(c1) = b(c1)-Flux(s)
        SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
              / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )
        A12 = DENs * SMDPN * volume(c1) / A % sav(c1)
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
!>>>>> write(*,*) 'mass before', errmax !, u(cm), w(cm), p(cm)

!------------------------------------------------!
!     Solve the pressure correction equation     !
!------------------------------------------------!
!>>>>> A % val(1) = A % val(1) * 2.0

!---- Give the "false" flux back and set it to zero ! 2mat
!  do s=1,NS                                         ! 2mat
!    c1=SideC(1,s)                                   ! 2mat
!    c2=SideC(2,s)                                   ! 2mat
!    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat
!      if(StateMat(material(c1))==SOLID .or. &       ! 2mat
!         StateMat(material(c2))==SOLID) then        ! 2mat
!        b(c1) = b(c1) + Flux(s)                     ! 2mat
!        if(c2 > 0) b(c2) = b(c2) - Flux(s)          ! 2mat
!        Flux(s) = 0.0                               ! 2mat
!      end if                                        ! 2mat
!    end if                                          ! 2mat
!  end do                                            ! 2mat

!---- Disconnect the SOLID cells from FLUID system  ! 2mat
!  do s=1,NS                                         ! 2mat
!    c1=SideC(1,s)                                   ! 2mat
!    c2=SideC(2,s)                                   ! 2mat
!    if(c2>0 .or. c2<0.and.TypeBC(c2)==BUFFER) then  ! 2mat
!      if(c2 > 0) then ! => not BUFFER               ! 2mat
!        if(StateMat(material(c1)) == SOLID) then    ! 2mat
!          A12 = -A % val(A % pos(2,s))                  ! 2mat
!          A % val(A % pos(1,s)) = 0.0                   ! 2mat
!          A % val(A % pos(2,s)) = 0.0                   ! 2mat
!          if(StateMat(material(c2)) == FLUID) then  ! 2mat
!            A % val(A % dia(c2)) = A % val(A % dia(c2)) -  A12  ! 2mat
!          endif                                     ! 2mat
!        end if                                      ! 2mat
!        if(StateMat(material(c2)) == SOLID) then    ! 2mat
!          A12 = -A % val(A % pos(1,s))                  ! 2mat
!          A % val(A % pos(2,s)) = 0.0                   ! 2mat
!          A % val(A % pos(1,s)) = 0.0                   ! 2mat
!          if(StateMat(material(c1)) == FLUID) then  ! 2mat
!            A % val(A % dia(c1)) = A % val(A % dia(c1)) -  A12  ! 2mat
!          endif                                     ! 2mat
!        end if                                      ! 2mat
!      else            ! => BUFFER                   ! 2mat
!        if(StateMat(material(c1)) == SOLID  .or. &  ! 2mat
!           StateMat(material(c2)) == SOLID) then    ! 2mat
!          A12 = -A % bou(c2)                           ! 2mat
!          A % bou(c2) = 0.0                            ! 2mat
!          if(StateMat(material(c1)) == FLUID) then  ! 2mat
!            A % val(A % dia(c1)) = A % val(A % dia(c1)) -  A12  ! 2mat
!          endif                                     ! 2mat
!        end if                                      ! 2mat
!      end if                                        ! 2mat
!    end if                                          ! 2mat
!  end do                                            ! 2mat

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!     Don't solve the pressure corection too accurate.      !
!     Value 1.e-18 blows the solution.                      !
!     Value 1.e-12 keeps the solution stable                !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!  
  niter=40
  call cgs(NC, NbC, A,                     &
        PP % n, b, PREC, niter, PP % STol,  &
        res(4), error) 
  write(LineRes(53:64),  '(1PE12.3)') res(4)
  write(LineRes(89:92),  '(I4)')      niter
!-----------------------------------!
!     Update the pressure field     !
!-----------------------------------!
  do c = 1, NC
    P % n(c)  =  P % n(c)  +  P % URF  *  PP % n(c)
  end do
!--------------------------------------!
!     Normalize the pressure field     !
!--------------------------------------!
!  Pmax  = maxval(P % n(1:NC))
!  Pmin  = minval(P % n(1:NC))
 
!  call glomax(Pmax)
!  call glomin(Pmin)
 
!bn003  P % n = P % n - 0.5*(Pmax+Pmin)
 
  call Exchng(PP % n)

  RETURN 

  end subroutine CalcPS
