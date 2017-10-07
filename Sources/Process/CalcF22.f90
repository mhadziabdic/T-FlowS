!======================================================================!
  subroutine CalcF22(var, phi,             &
                      dphidx, dphidy, dphidz)
!----------------------------------------------------------------------!
! Discretizes and solves eliptic relaxation equations for f22          !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer       :: var
  type(Unknown) :: phi
  real          :: dphidx(-NbC:NC), dphidy(-NbC:NC), dphidz(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  integer :: s, c, c1, c2, niter, miter
  real    :: Fex, Fim 
  real    :: A0, A12, A21
  real    :: error
  real    :: dphidxS, dphidyS, dphidzS
!======================================================================! 
!  The form of equations which are solved:
!
!     /            /              /
!    |  df22      | f22 dV       | f22hg dV
!  - |  ---- dS + | ------   =   | --------
!    |   dy       |  Lsc^2       |  Lsc^2
!   /            /              /
!
!
!  Dimension of the system under consideration
!
!     [A]{f22} = {b}   [kg K/s]
!
!  Dimensions of certain variables:
!     f22            [1/s]
!     Lsc            [m]
!
!======================================================================!

  Aval = 0.0

  b=0.0


!----- This is important for "copy" boundary conditions. Find out why !
  do c=-NbC,-1
    Abou(c)=0.0
  end do

!-----------------------------------------! 
!     Initialize variables and fluxes     !
!-----------------------------------------! 

!----- old values (o and oo)
  if(ini == 1) then
    do c=1,NC
      phi % oo(c)  = phi % o(c)
      phi % o (c)  = phi % n(c)
      phi % Doo(c) = phi % Do(c)
      phi % Do (c) = 0.0 
      phi % Xoo(c) = phi % Xo(c)
      phi % Xo (c) = phi % X(c) 
    end do
  end if


!----------------------------------------------------!
!     Browse through all the faces, where else ?     !
!----------------------------------------------------!

!----- new values
  do c=1,NC
    phi % X(c) = 0.0
  end do

!==================!
!                  !
!     Difusion     !
!                  !
!==================!

!--------------------------------!
!     Spatial Discretization     !
!--------------------------------!
  do s=1,NS       

    c1=SideC(1,s)
    c2=SideC(2,s)   

    dphidxS = fF(s)*dphidx(c1) + (1.0-fF(s))*dphidx(c2)
    dphidyS = fF(s)*dphidy(c1) + (1.0-fF(s))*dphidy(c2)
    dphidzS = fF(s)*dphidz(c1) + (1.0-fF(s))*dphidz(c2)


!---- total (exact) diffusive flux
    Fex=( dphidxS*Sx(s) + dphidyS*Sy(s) + dphidzS*Sz(s) )

    A0 =  Scoef(s)

!---- implicit diffusive flux
!.... this_proc is a very crude approximation: Scoef is not
!.... corrected at interface between materials
    Fim=( dphidxS*Dx(s)                      &
         +dphidyS*Dy(s)                      &
         +dphidzS*Dz(s))*A0

!---- this_proc is yet another crude approximation:
!.... A0 is calculated approximatelly
!    if( StateMat(material(c1))==FLUID .and.  &  ! 2mat
!        StateMat(material(c2))==SOLID        &  ! 2mat
!        .or.                                 &  ! 2mat 
!        StateMat(material(c1))==SOLID .and.  &  ! 2mat
!        StateMat(material(c2))==FLUID ) then    ! 2mat
!      A0 = A0 + A0                              ! 2mat
!    end if                                      ! 2mat

!---- straight diffusion part 
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

!---- cross diffusion part
    phi % X(c1) = phi % X(c1) + Fex - Fim 
    if(c2  > 0) then
      phi % X(c2) = phi % X(c2) - Fex + Fim 
    end if 

!----- calculate the coefficients for the sysytem matrix
    if( (DIFFUS == CN) .or. (DIFFUS == FI) ) then  

      if(DIFFUS  ==  CN) then       ! Crank Nicholson
        A12 = 0.5 * A0 
        A21 = 0.5 * A0 
      end if

      if(DIFFUS  ==  FI) then       ! Fully implicit
        A12 = A0 
        A21 = A0
      end if

!----- fill the system matrix
      if(c2  > 0) then
        Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12
        Aval(Adia(c1))    = Aval(Adia(c1))    + A12
        Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21
        Aval(Adia(c2))    = Aval(Adia(c2))    + A21
      else if(c2  < 0) then
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Outflow is not included because it was causing problems     !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -! 
        if( (TypeBC(c2) == INFLOW)) then                    
!---------  (TypeBC(c2) == OUTFLOW) ) then   
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          b(c1) = b(c1) + A12 * phi % n(c2)

        else

        if( (TypeBC(c2) == WALL).or.                          &
            (TypeBC(c2) == WALLFL) ) then
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
!=============================================================!
! Source coefficient is filled in SourceF22.f90 in order to   !
! get updated values of f22 on the wall.                      !
! Othrwise f22 equation does not converge very well           !
!          b(c1) = b(c1) + A12 * phi % n(c2)                  !
!=============================================================!
        else if( TypeBC(c2) == BUFFER ) then  
          Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          Abou(c2) = - A12  ! cool parallel stuff
        endif
      end if     
     end if
    end if

  end do  ! through sides

!---------------------------------!
!     Temporal discretization     !
!---------------------------------!

!----- Adams-Bashfort scheeme for diffusion fluxes
  if(DIFFUS == AB) then 
    do c=1,NC
      b(c) = b(c) + 1.5 * phi % Do(c) - 0.5 * phi % Doo(c)
    end do  
  end if

!----- Crank-Nicholson scheme for difusive terms
  if(DIFFUS == CN) then 
    do c=1,NC
      b(c) = b(c) + 0.5 * phi % Do(c)
    end do  
  end if
                 
!----- Fully implicit treatment for difusive terms
!      is handled via the linear system of equations 

!----- Adams-Bashfort scheeme for cross diffusion 
  if(CROSS == AB) then
    do c=1,NC
      b(c) = b(c) + 1.5 * phi % Xo(c) - 0.5 * phi % Xoo(c)
    end do 
  end if

!----- Crank-Nicholson scheme for cross difusive terms
  if(CROSS == CN) then
    do c=1,NC
      b(c) = b(c) + 0.5 * phi % X(c) + 0.5 * phi % Xo(c)
    end do 
  end if

!----- Fully implicit treatment for cross difusive terms
  if(CROSS == FI) then
    do c=1,NC
      b(c) = b(c) + phi % X(c)
    end do 
  end if

!========================================!
!                                        !  
!     Source terms and wall function     !
!     (Check if it is good to call it    !
!      before the under relaxation ?)    !
!                                        !
!========================================!

  if(SIMULA == EBM) then
    call SourceF22_EBM
  else
    call SourceF22KEPSV2F()
  end if

!=====================================!
!                                     !
!     Solve the equations for phi     !
!                                     !    
!=====================================!
    do c=1,NC
      b(c) = b(c) + Aval(Adia(c)) * (1.0-phi % URF)*phi % n(c) / phi % URF
      Aval(Adia(c)) = Aval(Adia(c)) / phi % URF
    end do 


  if(ALGOR == SIMPLE)   miter=300
  if(ALGOR == FRACT)    miter=5

  niter=miter
  call cg(NC, Nbc, NONZERO, Aval,Acol,Arow,Adia,Abou,   &
           phi % n, b, PREC,                            &
           niter,phi % STol, res(var), error)

  
  if(this_proc < 2) write(*,*) 'Var ', var, res(var), niter 

  call Exchng(phi % n)

  RETURN

  end subroutine CalcF22
