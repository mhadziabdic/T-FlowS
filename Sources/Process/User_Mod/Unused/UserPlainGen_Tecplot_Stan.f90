!======================================================================!
  subroutine UserPlainGen_Tecplot_Stan(grid, n)
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
!-----------------------------[Arguments]------------------------------!
    integer :: n
    real    :: x1_dis, x0_dis, y1_dis, y0_dis, z1_dis, z0_dis
    real    :: TauWup, TauWdown, r23, TKE, Cmu_mod, Pk_tmp, Cs_mod, Cs_prime
    real    :: Utan, UnorSq, Unor, UtotSq, dely, Stot, R, UtanX, TKE_mean
    real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32, a_mn_a_mn
    real    :: u11, u22, u33, u12, u13, u21, u31, u23, u32
    real    :: S11, S22, S33, S12, S13, S21, S31, S23, S32, a_lk_s_lk
!-------------------------------[Locals]-------------------------------!
    integer             :: c, s, c1, c2
    character*24        :: inflowfile
    character*80        :: namout
    character*39        :: path
!======================================================================!
   TauWdown = 0.0
   TauWup   = 0.0 
  
!  open(19, file='Cyl_data.dat', status='OLD', position = 'append')

  open(9, file='Slice.dat')
  if(this < 2) print *, '# Now reading the file: Slice.dat '
  read(9,*) x1_dis, x0_dis
  read(9,*) y1_dis, y0_dis
  read(9,*) z1_dis, z0_dis
  if(this < 2) print *, '# X:[ ', x0_dis, " ;", x1_dis, "]", &
  ' Y:[ ', y0_dis, " ;", y1_dis, "]", ' Z:[ ', z0_dis, " ;", z1_dis, "]"

    inflowfile = 'inflow-xxxxxx-xxx.dat'

    write(inflowfile(8:13),'(I6.6)') n
    write(inflowfile(15:17),'(I3.3)') this
    open(500+this, file=inflowfile)

    if( this < 2 ) then
      print *, 'Capturing field..'
    end if


!     call GradP(P % n,Px,Py,Pz)
!     call CalcShear(U%mean,V%mean,W%mean,Shear)
!     call CalcVort()


!    call GraPhi(U % mean, 1, Ux,.TRUE.)    ! dU/dx
!    call GraPhi(U % mean, 2, Uy,.TRUE.)    ! dU/dy
!    call GraPhi(U % mean, 3, Uz,.TRUE.)    ! dU/dz

!    call GraPhi(V % mean, 1, Vx,.TRUE.)    ! dV/dx
!    call GraPhi(V % mean, 2, Vy,.TRUE.)    ! dV/dy
!    call GraPhi(V % mean, 3, Vz,.TRUE.)    ! dV/dz

!    call GraPhi(W % mean, 1, Wx,.TRUE.)    ! dW/dx
!    call GraPhi(W % mean, 2, Wy,.TRUE.)    ! dW/dy
!    call GraPhi(W % mean, 3, Wz,.TRUE.)    ! dW/dz

    do c=1, NC
      if( xc(c)<x1_dis.and.xc(c)>x0_dis.and.yc(c)<y1_dis.and. &
          yc(c)>y0_dis.and.zc(c)<z1_dis.and.zc(c)>z0_dis) then
!        yc(c) > y0_dis .and. WallDs(c) < z1_dis .and. WallDs(c) > z0_dis) then
        if (SIMULA == ZETA) then
    r23 = 2.0/3.0    
    TKE = 0.5*(uu%mean(c)+vv%mean(c)+ww%mean(c)) + &
    0.5*(VAR10x(c)-U%mean(c)*U%mean(c)+VAR10y(c)-V%mean(c)*V%mean(c)+VAR10z(c)-W%mean(c)*W%mean(c))

    a11 = (uu%mean(c)+VAR10x(c)-U%mean(c)*U%mean(c))/TKE - r23 
    a22 = (vv%mean(c)+VAR10y(c)-V%mean(c)*V%mean(c))/TKE - r23 
    a33 = (ww%mean(c)+VAR10z(c)-W%mean(c)*W%mean(c))/TKE - r23 
    a12 = (uv%mean(c)+VAR11x(c)-U%mean(c)*V%mean(c))/TKE 
    a21 = a12 
    a13 = (uw%mean(c)+VAR11y(c)-U%mean(c)*W%mean(c))/TKE 
    a31 = a13
    a23 = (vw%mean(c)+VAR11z(c)-V%mean(c)*W%mean(c))/TKE 
    a32 = a23

   S11 = Ux(c)
   S22 = Vy(c)
   S33 = Wz(c)
   S12 = 0.5*(Uy(c)+Vx(c))
   S21 = S12
   S13 = 0.5*(Uz(c)+Wx(c))
   S31 = S13
   S23 = 0.5*(Vz(c)+Wy(c))
   S32 = S23

   a_mn_a_mn = sqrt(a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23))
   a_lk_s_lk = a11*S11 + a22*S22 + a33*S33 + 2.0*(a12*S12+a13*S13+a23*S23)

   Pk_tmp = -(a11*Ux(c) + a12*Uy(c) + a13*Uz(c) +&
              a12*Vx(c) + a22*Vy(c) + a23*Vz(c) +&
              a13*Wx(c) + a23*Wy(c) + a33*Wz(c))

!   Pk_tmp = Pk_tmp/(a_mn_a_mn*Shear(c))

           write (500+this,'(7E17.7E3)') xc(c), zc(c), U % n(c), V % n(c), W % n(c), Kin % n(c), Pk_tmp 
!           write (500+this,'(19E17.7E3)') xc(c), yc(c), U % mean(c), V % mean(c), W % mean(c), P % mean(c), &
!                                                  UU % mean(c) - U % mean(c)*U % mean(c), & !uu resolved
!                                                  VV % mean(c) - V % mean(c)*V % mean(c), & !vv resolved
!                                                  Ww % mean(c) - W % mean(c)*W % mean(c), & !ww resolved
!                                                  UV % mean(c) - U % mean(c)*V % mean(c), & !uv resolved
!                                                  UW % mean(c) - U % mean(c)*W % mean(c), & !uw resolved
!                                                  VW % mean(c) - V % mean(c)*W % mean(c), & !vw resolved
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Ux(c) ), & !uu modelled
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Vy(c) ), & !vv modelled
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Wz(c) ), & !ww modelled
!                                                  - VISt_mean(c)*( Uy(c) + Vx(c) ), & !uv modelled
!                                                  - VISt_mean(c)*( Uz(c) + Wx(c) ), & !uw modelled
!                                                  - VISt_mean(c)*( Vz(c) + Wy(c) ), & !vw modelled
!                                                  Eps % mean(c)
        elseif (SIMULA == HJ) then

    r23 = 2.0/3.0    
    TKE = 0.5*(uu%n(c)+vv%n(c)+ww%n(c)) 

    a11 = uu%n(c)/TKE - r23 
    a22 = vv%n(c)/TKE - r23 
    a33 = ww%n(c)/TKE - r23 
    a12 = uv%n(c)/TKE 
    a21 = a12 
    a13 = uw%n(c)/TKE 
    a31 = a13
    a23 = vw%n(c)/TKE 
    a32 = a23
    u11 = uu%n(c)
    u22 = vv%n(c)
    u33 = ww%n(c)
    u12 = uv%n(c)
    u21 = a12 
    u13 = uw%n(c)
    u31 = a13
    u23 = vw%n(c)
    u32 = a23

   S11 = Ux(c)
   S22 = Vy(c)
   S33 = Wz(c)
   S12 = 0.5*(Uy(c)+Vx(c))
   S21 = S12
   S13 = 0.5*(Uz(c)+Wx(c))
   S31 = S13
   S23 = 0.5*(Vz(c)+Wy(c))
   S32 = S23

   a_mn_a_mn = sqrt(a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23))
   a_lk_s_lk = a11*S11 + a22*S22 + a33*S33 + 2.0*(a12*S12+a13*S13+a23*S23)

   Pk_tmp = -(a11*Ux(c) + a12*Uy(c) + a13*Uz(c) +&
              a12*Vx(c) + a22*Vy(c) + a23*Vz(c) +&
              a13*Wx(c) + a23*Wy(c) + a33*Wz(c))

   Pk_tmp = -a_lk_s_lk
!   Pk_tmp2 =-(u11*Ux(c) + u12*Uy(c) + u13*Uz(c) +&
!              u12*Vx(c) + u22*Vy(c) + u23*Vz(c) +&
!              u13*Wx(c) + u23*Wy(c) + u33*Wz(c))

   Cs_prime = Pk_tmp/(a_mn_a_mn*Shear(c))
   Cs_mod   = Pk_tmp/(Shear(c))

   Cmu_mod = max(-(uu%n(c)*Ux(c)+vv%n(c)*Vy(c)+ww%n(c)*Wz(c)+&
                   uv%n(c)*(Vx(c)+Uy(c))+uw%n(c)*(Uz(c)+&
   Wx(c))+vw%n(c)*(Vz(c)+Wy(c)))/max(Kin%n(c)*Kin%n(c)/Eps%n(c)*Shear(c)*Shear(c),1.0e-12),0.0)

    write (500+this,'(9E17.7E3)') xc(c), zc(c), U % n(c), V % n(c), W%n(c), Cmu_mod, a_mn_a_mn, Cs_mod, Cs_prime 
!           write (500+this,'(30E17.7E3)') xc(c), yc(c), U % n(c), V % n(c), W%n(c), P % n(c), Kin%n(c),&
!                                          Kin%mean(c),&  
!                                        (VAR10x(c)-U%mean(c)*U%mean(c)),&
!                                        (VAR10y(c)-V%mean(c)*V%mean(c)),&
!                                        (VAR10z(c)-W%mean(c)*W%mean(c)),&    
!                                        (VAR11x(c)-U%mean(c)*V%mean(c)),&
!                                        (VAR11y(c)-U%mean(c)*W%mean(c)),&
!                                        (VAR11z(c)-V%mean(c)*W%mean(c)),&
!                                        uu % mean(c), & !uu modelled
!                                        vv % mean(c), & !vv modelled
!                                        ww % mean(c), & !ww modelled
!!                                        uv % mean(c), &!uv modelled
!                                        uw % mean(c), & !uw modelled
!                                        vw % mean(c),&
!                                        uu % n(c), & !uu modelled
!                                        vv % n(c), & !vv modelled
!                                        ww % n(c), & !ww modelled
!                                        uv % n(c), & !uv modelled
!!!!                                        uw % n(c), & !uw modelled
!                                        vw % n(c), &
!                                        (0.5*((VAR10x(c)-U%mean(c)*U%mean(c))+&
!                                             (VAR10y(c)-V%mean(c)*V%mean(c))+&
!                                             (VAR10z(c)-W%mean(c)*W%mean(c))))/TKE,&
!                                        (0.5*(uu%mean(c)+vv%mean(c)+ww%mean(c)))/TKE, &
!                                        Cmu_mod, VAR8x(c)!0.5*(Shear(c)-Vort(c))

        elseif(SIMULA == K_EPS) then
           write (500+this,'(10E17.7E3)') xc(c), yc(c), zc(c), U%n(c), V%n(c), W%n(c), P%n(c), Kin%n(c),&
                                          Eps%n(c), sqrt(Kin%n(c))*Cmu25*WallDs(c)/VISc  
        elseif (SIMULA == LES) then
           write (500+this,'(6E17.7E3)') xc(c), zc(c), U % n(c), V % n(c), W%n(c), P % n(c)!, &
!           write (500+this,'(11E17.7E3)') xc(c), yc(c), U % n(c), V % n(c), W%n(c), P % n(c), &
!                                        (uu % mean(c) - U%mean(c)*U%mean(c)),&
!                                        (vv % mean(c) - V%mean(c)*V%mean(c)),&
!                                        (ww % mean(c) - W%mean(c)*W%mean(c)),&
!                                        (uv % mean(c) - U%mean(c)*V%mean(c)),&
!                                        0.5*(Shear(c)-Vort(c))
        end if
      endif
    enddo



!    do s=1,NS
!      c1=SideC(1,s)
!      c2=SideC(2,s)

!      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
!        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!
!          UtotSq = U % n(c1) * U % n(c1) &
!                 + V % n(c1) * V % n(c1) &
!                 + W % n(c1) * W % n(c1)
!          Unor = ( U % n(c1) * Sx(s)     &
!                 + V % n(c1) * Sy(s)    &
!                 + W % n(c1) * Sz(s) )   &
!               / sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))
!          UnorSq = Unor*Unor
!          if( UtotSq   >  UnorSq) then
!            Utan = sqrt(UtotSq - UnorSq)
!          else
!            Utan = TINY
!          end if

!          if(yc(c2)>0.0.and.xc(c2)<0.0) then
!            TauWup = TauWup + VISc*Utan/WallDs(c1)*sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + &
!            Sz(s)*Sz(s))
!          else if(yc(c2)<0.0.and.xc(c2)<0.0) then
!            TauWdown = TauWdown + VISc*Utan/WallDs(c1)*sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + &
!            Sz(s)*Sz(s))
!          end if         
!        end if
!      end if
!    end do  

!    call GloSum(TauWup)
!    call GloSum(TauWdown)
    call wait
     
    close(9)
    close(500+this)
    call wait
!    if ( this < 2 ) then
!      print *, 'It is done'
!      write(19,*) n, TauWdown, TauWup 
!    end if

  end subroutine UserPlainGen_Tecplot_Stan
