!======================================================================!
  subroutine SourcesHJ(var)
!----------------------------------------------------------------------!
! Purpose:                                                             !
! Calculate source terms for transport equations for Re stresses and   !
! dissipation for HJ model.                                            !  
! Authors: Muhamed Hadziabdic                                          !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use les_mod

!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, s, c1, c2, i, var,icont
  real    :: Prod, Diss, VAR_hom, VAR_wall, mag, VAR_tot, Esor
  real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32
  real    :: n1,n2,n3,AA2,AA3,AA,Ret,ff2,fd,FF1,CC,C1W,C2W,fw,uu_nn
  real    :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  real    :: Eps11,Eps12,Eps13,Eps21,Eps22,Eps23,Eps31,Eps32,Eps33
  real    :: VAR_211,VAR_212,VAR_213,VAR_222,VAR_223,VAR_233,phi2_nn
  real    :: VAR_hm1,VAR_wal1,VAR_wal2, Feps, CC1P
  real    :: fss,e11w,e12w,e13w,e22w,e23w,e33w,E2,E3,EE,CC1,CC2
  real    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uyx, Uzx, Uzy
  real    :: Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
  real    :: Wxx, Wyy, Wzz, Wxy, Wxz, Wyz
  real    :: Diss_wall, Diss_hom, r13, r23, T1, T2
  real    :: S11, S22, S33, S12, S13, S21, S31, S23, S32
  real    :: V11, V22, V33, V12, V13, V21, V31, V23, V32
  real    :: CC3, CC4, CC5, a_lk_s_lk, PHI_hom, a_mn_a_mn
  real    :: VAR1w_11, VAR1w_22, VAR1w_33, VAR1w_12, VAR1w_13, VAR1w_23
  real    :: VAR2w_11, VAR2w_22, VAR2w_33, VAR2w_12, VAR2w_13, VAR2w_23
  real    :: VAR1_11, VAR1_22, VAR1_33, VAR1_12, VAR1_13, VAR1_23
  real    :: VAR2_11, VAR2_22, VAR2_33, VAR2_12, VAR2_13, VAR2_23
  real    :: VARip2_11,VARip2_22,VARip2_33,VARip2_12,VARip2_13,VARip2_23
  real    :: P11, P22, P33, P12, P13, P23, Eps1, Eps2
  real    :: duu_dx,duu_dy,duu_dz,dvv_dx,dvv_dy,dvv_dz,dww_dx,dww_dy,dww_dz
  real    :: duv_dx,duv_dy,duv_dz,duw_dx,duw_dy,duw_dz,dvw_dx,dvw_dy,dvw_dz
  real    :: dUdx, dUdy, dUdz, dVdx, dVdy, dVdz, dWdx, dWdy, dWdz 
  real,allocatable :: Diss1(:)
!======================================================================!

  allocate (Diss1(1:NC)); Diss1 = 0.0
  EE = 0.5
  AA = 0.5

  do c=1,NC
    Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), 1.0e-7)
    Lsc(c)=  (Kin % n(c))**1.5/Eps % n(c)
    Tsc(c)=  Kin % n(c)/Eps % n(c)
  end do

  call GraPhi(Kin%n,1,VAR6x,.TRUE.)             ! dK/dx
  call GraPhi(Kin%n,2,VAR6y,.TRUE.)             ! dK/dy
  call GraPhi(Kin%n,3,VAR6z,.TRUE.)             ! dK/dz

  call GraPhi(VAR6x,1,VAR7x,.TRUE.)             ! dK/dx
  call GraPhi(VAR6y,2,VAR7y,.TRUE.)             ! dK/dy
  call GraPhi(VAR6z,3,VAR7z,.TRUE.)             ! dK/dz

  do c = 1, NC
    Eps_tot(c) = Eps%n(c) + 0.5*VISc*(VAR7x(c)+VAR7y(c)+VAR7z(c))
  end do

!====================================================================!  
! Below is one of versions of HJ model that required much more memory!
!====================================================================!  
  if(var == 23) then
    call GraPhi(uu%n,1,VAR3x, .TRUE.) ! duu/dx  
    call GraPhi(uu%n,2,VAR3y, .TRUE.) ! duu/dy  
    call GraPhi(uu%n,3,VAR3z, .TRUE.) ! duu/dz  

    call GraPhi(vv%n,1,VAR4x, .TRUE.) ! duw/dx  
    call GraPhi(vv%n,2,VAR4y, .TRUE.) ! duw/dy  
    call GraPhi(vv%n,3,VAR4z, .TRUE.) ! duw/dz  

    call GraPhi(ww%n,1,VAR5x, .TRUE.) ! duw/dx  
    call GraPhi(ww%n,2,VAR5y, .TRUE.) ! duw/dy  
    call GraPhi(ww%n,3,VAR5z, .TRUE.) ! duw/dz  

    call GraPhi(uv%n,1,VAR6x, .TRUE.) ! duv/dx  
    call GraPhi(uv%n,2,VAR6y, .TRUE.) ! duv/dy  
    call GraPhi(uv%n,3,VAR6z, .TRUE.) ! duv/dz  

    call GraPhi(uw%n,1,VAR7x, .TRUE.) ! duw/dx  
    call GraPhi(uw%n,2,VAR7y, .TRUE.) ! duw/dy  
    call GraPhi(uw%n,3,VAR7z, .TRUE.) ! duw/dz  

    call GraPhi(vw%n,1,VAR8x, .TRUE.) ! duw/dx  
    call GraPhi(vw%n,2,VAR8y, .TRUE.) ! duw/dy  
    call GraPhi(vw%n,3,VAR8z, .TRUE.) ! duw/dz  

    call GraPhi(Ux,1,VAR1x, .TRUE.)  ! d2U/dxdx
    call GraPhi(Uy,2,VAR1y, .TRUE.)  ! d2U/dydy
    call GraPhi(Uz,3,VAR1z, .TRUE.)  ! d2U/dzdz
    call GraPhi(Ux,2,VAR2x, .TRUE.)  ! d2U/dxdy
    call GraPhi(Ux,3,VAR2y, .TRUE.)  ! d2U/dxdz
    call GraPhi(Uy,3,VAR2z, .TRUE.)  ! d2U/dydz

    call GraPhi(Vx,1,VAR9x, .TRUE.)  ! d2V/dxdx
    call GraPhi(Vy,2,VAR9y, .TRUE.)  ! d2V/dydy
    call GraPhi(Vz,3,VAR9z, .TRUE.)  ! d2V/dzdz
    call GraPhi(Vx,2,VAR10x, .TRUE.)  ! d2V/dxdy
    call GraPhi(Vx,3,VAR10y, .TRUE.)  ! d2V/dxdz
    call GraPhi(Vy,3,VAR10z, .TRUE.)  ! d2V/dydz

    call GraPhi(Wx,1,VAR11x, .TRUE.)  ! d2W/dxdx
    call GraPhi(Wy,2,VAR11y, .TRUE.)  ! d2W/dydy
    call GraPhi(Wz,3,VAR11z, .TRUE.)  ! d2W/dzdz
    call GraPhi(Wx,2,VAR12x, .TRUE.)  ! d2W/dxdy
    call GraPhi(Wx,3,VAR12y, .TRUE.)  ! d2W/dxdz
    call GraPhi(Wy,3,VAR12z, .TRUE.)  ! d2W/dydz

    do c=1,NC
      Uxx = VAR1x(c)
      Uyy = VAR1y(c)
      Uzz = VAR1z(c)
      Uxy = VAR2x(c)
      Uxz = VAR2y(c)
      Uyz = VAR2z(c)
      Vxx = VAR9x(c)
      Vyy = VAR9y(c)
      Vzz = VAR9z(c)
      Vxy = VAR10x(c)
      Vxz = VAR10y(c)
      Vyz = VAR10z(c)
      Wxx = VAR11x(c)
      Wyy = VAR11y(c)
      Wzz = VAR11z(c)
      Wxy = VAR12x(c)
      Wxz = VAR12y(c)
      Wyz = VAR12z(c)
      dUdx= Ux(c) 
      dUdy= Uy(c) 
      dUdz= Uz(c) 
      dVdx= Vx(c) 
      dVdy= Vy(c) 
      dVdz= Vz(c) 
      dWdx= Wx(c) 
      dWdy= Wy(c) 
      dWdz= Wz(c) 
      duu_dx = VAR3x(c)  
      duu_dy = VAR3y(c)  
      duu_dz = VAR3z(c)  
      dvv_dx = VAR4x(c)  
      dvv_dy = VAR4y(c)  
      dvv_dz = VAR4z(c)  
      dww_dx = VAR5x(c)  
      dww_dy = VAR5y(c)  
      dww_dz = VAR5z(c)  
      duv_dx = VAR6x(c)  
      duv_dy = VAR6y(c)  
      duv_dz = VAR6z(c)  
      duw_dx = VAR7x(c)  
      duw_dy = VAR7y(c)  
      duw_dz = VAR7z(c)  
      dvw_dx = VAR8x(c)  
      dvw_dy = VAR8y(c)  
      dvw_dz = VAR8z(c)  

      Diss1(c) =  -2.0*VISc*(duu_dx*Uxx + duv_dy*Uyy + duw_dz*Uzz + &
                            Uxy*(duv_dx+duu_dy) + Uxz*(duw_dx+duu_dz) + Uyz*(duw_dy + duv_dz) +&
                            duv_dx*Vxx + dvv_dy*Vyy + dvw_dz*Vzz + &
                            Vxy*(dvv_dx+duv_dy) + Vxz*(dvw_dx+duv_dz) + Vyz*(dvw_dy + dvv_dz) +&
                            duw_dx*Wxx + dvw_dy*Wyy + dww_dz*Wzz + &
                            Wxy*(dvw_dx+duw_dy) + Wxz*(dww_dx+duw_dy) + Wyz*(dww_dy + dvw_dz) +&
    0.32*Kin%n(c)/Eps%n(c)*(Uxx    * (duu_dx*dUdx  + duv_dx*dUdy  + duw_dx*dUdz ) + & 
                            Uyy    * (duv_dy*dUdx  + dvv_dy*dUdy  + dvw_dy*dUdz ) + & 
                            Uzz    * (duw_dz*dUdx  + dvw_dz*dUdy  + dww_dz*dUdz ) + & 
                            Uxy    * (duu_dy*dUdx  + duv_dy*dUdy  + duw_dy*dUdz   + &
                                      duv_dx*dUdx  + dvv_dx*dUdy  + dvw_dx*dUdz ) + & 
                            Uxz    * (duu_dz*dUdx  + duv_dz*dUdy  + duw_dz*dUdz   + &
                                      duw_dx*dUdx  + dvw_dx*dUdy  + dww_dx*dUdz ) + &
                            Uyz    * (duv_dz*dUdx  + dvv_dz*dUdy  + dvw_dz*dUdz   + &
                                      duw_dy*dUdx  + dvw_dy*dUdy  + dww_dy*dUdz ) + &
                            Vxx    * (duu_dx*dVdx  + duv_dx*dVdy  + duw_dx*dVdz ) + & 
                            Vyy    * (duv_dy*dVdx  + dvv_dy*dVdy  + dvw_dy*dVdz ) + & 
                            Vzz    * (duw_dz*dVdx  + dvw_dz*dVdy  + dww_dz*dVdz ) + & 
                            Vxy    * (duu_dy*dVdx  + duv_dy*dVdy  + duw_dy*dVdz   + &
                                      duv_dx*dVdx  + dvv_dx*dVdy  + dvw_dx*dVdz ) + & 
                            Vxz    * (duu_dz*dVdx  + duv_dz*dVdy  + duw_dz*dVdz   + &
                                      duw_dx*dVdx  + dvw_dx*dVdy  + dww_dx*dVdz ) + &
                            Vyz    * (duv_dz*dVdx  + dvv_dz*dVdy  + dvw_dz*dVdz   + &
                                      duw_dy*dVdx  + dvw_dy*dVdy  + dww_dy*dVdz ) + &
                            Wxx    * (duu_dx*dWdx  + duv_dx*dWdy  + duw_dx*dWdz ) + & 
                            Wyy    * (duv_dy*dWdx  + dvv_dy*dWdy  + dvw_dy*dWdz ) + & 
                            Wzz    * (duw_dz*dWdx  + dvw_dz*dWdy  + dww_dz*dWdz ) + & 
                            Wxy    * (duu_dy*dWdx  + duv_dy*dWdy  + duw_dy*dWdz   + &
                                      duv_dx*dWdx  + dvv_dx*dWdy  + dvw_dx*dWdz ) + & 
                            Wxz    * (duu_dz*dWdx  + duv_dz*dWdy  + duw_dz*dWdz   + &
                                      duw_dx*dWdx  + dvw_dx*dWdy  + dww_dx*dWdz ) + &
                            Wyz    * (duv_dz*dWdx  + dvv_dz*dWdy  + dvw_dz*dWdz   + &
                                      duw_dy*dWdx  + dvw_dy*dWdy  + dww_dy*dWdz ))) 
    end do
  end if

  if(var == 13) then
  do i=1,3
    if(i == 1) then
      call GraPhi(Ux,1,VAR2x, .TRUE.)  ! d2U/dxdx
      call GraPhi(Ux,2,VAR2y, .TRUE.)  ! d2U/dxdy
      call GraPhi(Ux,3,VAR2z, .TRUE.)  ! d2U/dxdz
      call GraPhi(Uy,2,VAR1x, .TRUE.) ! d2U/dydy
      call GraPhi(Uy,3,VAR1y, .TRUE.) ! d2U/dydz
      call GraPhi(Uz,3,VAR1z, .TRUE.) ! d2U/dzdz
    end if
    if(i == 2) then
      call GraPhi(Vx,1,VAR2x, .TRUE.)  ! d2V/dxdx
      call GraPhi(Vx,2,VAR2y, .TRUE.)  ! d2V/dxdy
      call GraPhi(Vx,3,VAR2z, .TRUE.)  ! d2V/dxdz
      call GraPhi(Vy,2,VAR1x, .TRUE.) ! d2V/dydy
      call GraPhi(Vy,3,VAR1y, .TRUE.) ! d2V/dydz
      call GraPhi(Vz,3,VAR1z, .TRUE.) ! d2V/dzdz
    end if
    if(i == 3) then
      call GraPhi(Wx,1,VAR2x, .TRUE.)  ! d2W/dxdx
      call GraPhi(Wx,2,VAR2y, .TRUE.)  ! d2W/dxdy
      call GraPhi(Wx,3,VAR2z, .TRUE.)  ! d2W/dxdz
      call GraPhi(Wy,2,VAR1x, .TRUE.) ! d2W/dydy
      call GraPhi(Wy,3,VAR1y, .TRUE.) ! d2W/dydz
      call GraPhi(Wz,3,VAR1z, .TRUE.) ! d2W/dzdz
    end if

    do c=1,NC
      if(i == 1) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = &
                2.0*0.25*VISc*Kin%n(c)/Eps_tot(c)*&
               (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
      if(i == 2) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = Diss1(c) +&
                2.0*0.25*VISc*Kin%n(c)/Eps_tot(c)*&
                (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
      if(i == 3) then
        Uxx = VAR2x(c)
        Uxy = VAR2y(c)
        Uyx = Uxy
        Uxz = VAR2z(c)
        Uzx = Uxz
        Uyy = VAR1x(c)
        Uyz = VAR1y(c)
        Uzy = Uyz
        Uzz = VAR1z(c)
        Diss1(c) = Diss1(c) +&
                2.0*0.25*VISc*Kin%n(c)/Eps_tot(c)*&
                (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+&
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+&
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+&
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+&
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+&
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+&
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+&
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+&
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
    end do
  end do  ! i
  end if                               

  call GraPhi(Lsc,1,VAR2x,.TRUE.)             ! df22/dx
  call GraPhi(Lsc,2,VAR2y,.TRUE.)             ! df22/dy
  call GraPhi(Lsc,3,VAR2z,.TRUE.)             ! df22/dz

  r23 = 2.0/3.0 
  r13 = 1.0/3.0 
  do  c = 1, NC
    Pk(c)      = max(-(uu % n(c)*Ux(c) + uv % n(c)*Uy(c) + uw % n(c)*Uz(c) +&
                       uv % n(c)*Vx(c) + vv % n(c)*Vy(c) + vw % n(c)*Vz(c) +&
                       uw % n(c)*Wx(c) + vw % n(c)*Wy(c) + ww % n(c)*Wz(c)),1.0e-10)                
  
    mag = max(0.0,sqrt(VAR2x(c)*VAR2x(c)+VAR2y(c)*VAR2y(c)+VAR2z(c)*VAR2z(c)),1.0e-10)       
    n1 = VAR2x(c)/mag 
    n2 = VAR2y(c)/mag
    n3 = VAR2z(c)/mag

    a11 = uu % n(c)/Kin % n(c) - r23 
    a22 = vv % n(c)/Kin % n(c) - r23
    a33 = ww % n(c)/Kin % n(c) - r23
    a12 = uv % n(c)/Kin % n(c)   
    a21 = a12
    a13 = uw % n(c)/Kin % n(c)    
    a31 = a13
    a23 = vw % n(c)/Kin % n(c)    
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

    V11 = 0.0
    V22 = 0.0
    V33 = 0.0
    V12 = 0.5*(Uy(c)-Vx(c)) - omegaZ
    V21 = -V12 + omegaZ
    V13 = 0.5*(Uz(c)-Wx(c)) + omegaY
    V31 = -V13 - omegaY
    V23 = 0.5*(Vz(c)-Wy(c)) - omegaX
    V32 = -V23 + omegaX

    AA2=(a11**2)+(a22**2)+(a33**2)+2*((a12**2)+(a13**2)+(a23**2))

    AA3= a11**3 + a22**3 + a33**3 + &
         3*a12**2*(a11+a22) + 3*a13**2*(a11+a33) +&
         3*a23**2*(a22+a33) + 6*a12*a13*a23

    AA=1.0 - (9.0/8.0)*(AA2-AA3)
    AA=max(AA,0.0)
    AA=min(AA,1.0)
 
    uu_nn     = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3 &
               + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3 &
               + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)

    a_mn_a_mn = a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23)
    a_lk_s_lk = a11*S11 + a22*S22 + a33*S33 + 2.0*(a12*S12+a13*S13+a23*S23)
 
    EE=AA
    fss=1.0-(AA**0.5*EE**2.0)
    do icont=1,6
      Eps11= (1.0 - fss)*r23*Eps %n(c) + fss * uu%n(c)/Kin%n(c) * Eps%n(c)
      Eps22= (1.0 - fss)*r23*Eps %n(c) + fss * vv%n(c)/Kin%n(c) * Eps%n(c)
      Eps33= (1.0 - fss)*r23*Eps %n(c) + fss * ww%n(c)/Kin%n(c) * Eps%n(c)
      Eps12= fss * uv%n(c)/Kin%n(c) * Eps%n(c)
      Eps13= fss * uw%n(c)/Kin%n(c) * Eps%n(c)
      Eps23= fss * vw%n(c)/Kin%n(c) * Eps%n(c) 
      Eps21= Eps12
      Eps31= Eps13
      Eps32= Eps23

      e11= Eps11/Eps%n(c) - r23
      e22= Eps22/Eps%n(c) - r23
      e33= Eps33/Eps%n(c) - r23
      e12= Eps12/Eps%n(c)
      e13= Eps13/Eps%n(c)
      e23= Eps23/Eps%n(c)
      e21= e12
      e31= e13
      e32= e23
      E2=(e11**2)+(e22**2)+(e33**2)+2*((e12**2)+(e13**2)+(e23**2))

      E3= e11**3 + e22**3 + e33**3 + &
         3*e12**2*(e11+e22) + 3*e13**2*(e11+e33) +&
         3*e23**2*(e22+e33) + 6*e12*e13*e23

      EE=1.0 - (9.0/8.0)*(E2-E3)

      EE=max(EE,0.0)
      EE=min(EE,1.0)
      fss=1.0-(AA**0.5*EE**2.0)
    end do
     
    Ret= (Kin % n(c)*Kin % n(c))/(VISc*Eps_tot(c)+tiny)
    Feps = 1.0 - ((Ce2-1.4)/Ce2)*exp(-(Ret/6.0)**2.0)
    ff2=min((Ret/150)**1.5, 1.0)
    fd=1.0/(1.0+0.1*Ret)
    FF1=min(0.6, AA2)
    CC=2.5*AA*FF1**0.25*ff2
    CC1=CC+SQRT(AA)*(EE**2)
    CC2=0.8*SQRT(AA)
    C1W=max((1.0-0.7*CC), 0.3)
    C2W=min(AA,0.3)
    fw=min((Kin%n(c)**1.5)/(2.5*Eps_tot(c)*WallDs(c)),1.4)


    P11     = -2.0*(uu %n(c)*Ux(c) + uv % n(c)*Uy(c) + uw % n(c)*Uz(c)) &
              -2.0*omegaY*2.0*uw%n(c) + 2.0*omegaZ*2.0*uv%n(c) 
    P22     = -2.0*(uv %n(c)*Vx(c) + vv % n(c)*Vy(c) + vw % n(c)*Vz(c)) &
              +2.0*omegaX*2.0*vw%n(c) - 2.0*omegaZ*2.0*uw%n(c) 
    P33     = -2.0*(uw %n(c)*Wx(c) + vw % n(c)*Wy(c) + ww % n(c)*Wz(c)) &
              -2.0*omegaX*2.0*vw%n(c) + 2.0*omegaY*2.0*uw%n(c) 
    P12     = -(uu % n(c)*Vx(c) + uv % n(c)*Vy(c) + uw % n(c)*Vz(c) + &
                uv % n(c)*Ux(c) + vv % n(c)*Uy(c) + vw % n(c)*Uz(c)) &
                +2.0*omegaX*uw%n(c)-2.0*omegaY*vw%n(c)+2.0*omegaZ*(vv%n(c)-uu%n(c)) 
    P13     = -(uw % n(c)*Ux(c) + vw % n(c)*Uy(c) + ww % n(c)*Uz(c) + &
                uu % n(c)*Wx(c) + uv % n(c)*Wy(c) + uw % n(c)*Wz(c)) &
                -2.0*omegaX*uv%n(c)-2.0*omegaY*(ww%n(c)-uu%n(c))+2.0*omegaZ*vw%n(c) 
    P23     = -(uv % n(c)*Wx(c) + vv % n(c)*Wy(c) + vw % n(c)*Wz(c) + &
                uw % n(c)*Vx(c) + vw % n(c)*Vy(c) + ww % n(c)*Vz(c)) &
                -2.0*omegaX*(vw%n(c)-ww%n(c))+2.0*omegaY*uv%n(c)-2.0*omegaZ*uw%n(c) 

    VAR1_11 = -CC1*Eps%n(c)*a11 
    VAR1_22 = -CC1*Eps%n(c)*a22 
    VAR1_33 = -CC1*Eps%n(c)*a33 
    VAR1_12 = -CC1*Eps%n(c)*a12 
    VAR1_13 = -CC1*Eps%n(c)*a13 
    VAR1_23 = -CC1*Eps%n(c)*a23 
                                

    VAR2_11 = -CC2*(P11 - r23*Pk(c))
    VAR2_22 = -CC2*(P22 - r23*Pk(c))
    VAR2_33 = -CC2*(P33 - r23*Pk(c))
    VAR2_12 = -CC2*P12
    VAR2_13 = -CC2*P13
    VAR2_23 = -CC2*P23

    phi2_nn = VAR2_11*n1*n1+2*VAR2_12*n1*n2+2*VAR2_13*n1*n3+VAR2_22*n2*n2+2*VAR2_23*n2*n3+VAR2_33*n3*n3  

    VAR1w_11 = C1W*fw*Eps%n(c)/Kin%n(c)*(uu_nn-1.5*2.0*(uu%n(c)*n1*n1*0.0+uv%n(c)*n1*n2+uw%n(c)*n1*n3))
    VAR1w_22 = C1W*fw*Eps%n(c)/Kin%n(c)*(uu_nn-1.5*2.0*(uv%n(c)*n2*n1+vv%n(c)*n2*n2*0.0+vw%n(c)*n2*n3))
    VAR1w_33 = C1W*fw*Eps%n(c)/Kin%n(c)*(uu_nn-1.5*2.0*(uw%n(c)*n3*n1+vw%n(c)*n3*n2+ww%n(c)*n3*n3*0.0))
    VAR1w_12 = C1W*fw*Eps%n(c)/Kin%n(c)*(-1.5*(uu%n(c)*n2*n1+uv%n(c)*n2*n2*0.0+uw%n(c)*n2*n3 +&
                                               uv%n(c)*n1*n1*0.0+vv%n(c)*n1*n2+vw%n(c)*n1*n3)) 
    VAR1w_13 = C1W*fw*Eps%n(c)/Kin%n(c)*(-1.5*(uu%n(c)*n3*n1+uv%n(c)*n3*n2+uw%n(c)*n3*n3*0.0 +&
                                               uw%n(c)*n1*n1*0.0+vw%n(c)*n1*n2+ww%n(c)*n1*n3))
    VAR1w_23 = C1W*fw*Eps%n(c)/Kin%n(c)*(-1.5*(uw%n(c)*n2*n1+vw%n(c)*n2*n2*0.0+ww%n(c)*n2*n3 +&
                                               uv%n(c)*n3*n1+vv%n(c)*n3*n2+vw%n(c)*n3*n3)*0.0)

    VAR2w_11 = C2W*fw*(phi2_nn-1.5*2.0*(VAR2_11*n1*n1+VAR2_12*n1*n2+VAR2_13*n1*n3))
    VAR2w_22 = C2W*fw*(phi2_nn-1.5*2.0*(VAR2_12*n1*n2+VAR2_22*n2*n2+VAR2_23*n3*n2))
    VAR2w_33 = C2W*fw*(phi2_nn-1.5*2.0*(VAR2_13*n1*n3+VAR2_23*n2*n3+VAR2_33*n3*n3))
    VAR2w_12 = C2W*fw*(-1.5*(VAR2_11*n2*n1+VAR2_12*n2*n2+VAR2_13*n2*n3 +&
                             VAR2_12*n1*n1+VAR2_22*n1*n2+VAR2_23*n1*n3))
    VAR2w_13 = C2W*fw*(-1.5*(VAR2_11*n3*n1+VAR2_12*n3*n2+VAR2_13*n3*n3 +&
                             VAR2_13*n1*n1+VAR2_23*n1*n2+VAR2_33*n1*n3))
    VAR2w_23 = C2W*fw*(-1.5*(VAR2_13*n2*n1+VAR2_23*n2*n2+VAR2_33*n2*n3 +&
                             VAR2_12*n3*n1+VAR2_22*n3*n2+VAR2_23*n3*n3))

!------- uu stress
    if(var == 6) then
!==============================================================================================================================!
      b(c) = b(c) + (max(P11,0.0)+CC1*Eps%n(c)*r23+max(VAR2_11,0.0)+max(VAR1w_11,0.0)+max(VAR2w_11,0.0))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*3.0*n1*n1 + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P11,0.0)+max(-VAR2_11,0.0)+max(-VAR1w_11,0.0)+max(-VAR2w_11,0.0) + &
                      (1.0-fss)*r23*Eps%n(c))/max(uu%n(c),1.0e-10)*volume(c) 
!==============================================================================================================================!
!------- vv stress
    else if(var == 7) then
!==============================================================================================================================!
      b(c) = b(c) + (max(P22,0.0)+CC1*Eps%n(c)*r23+max(VAR2_22,0.0)+max(VAR1w_22,0.0)+max(VAR2w_22,0.0))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*3.0*n2*n2 + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P22,0.0)+max(-VAR2_22,0.0)+max(-VAR1w_22,0.0)+max(-VAR2w_22,0.0)+ &
                      (1.0-fss)*r23*Eps%n(c))/max(vv%n(c),1.0e-10)*volume(c) 
!==============================================================================================================================!
!------- ww stress
    else if(var == 8) then
!==============================================================================================================================!
      b(c) = b(c) + (max(P33,0.0)+CC1*Eps%n(c)*r23+max(VAR2_33,0.0)+max(VAR1w_33,0.0)+max(VAR2w_33,0.0))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*3.0*n3*n3 + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P33,0.0)+max(-VAR2_33,0.0)+max(-VAR1w_33,0.0)+max(-VAR2w_33,0.0)+ &
                      (1.0-fss)*r23*Eps%n(c))/max(ww%n(c),1.0e-10)*volume(c) 
!==============================================================================================================================!
!==============================================================================================================================!
!------- uv stress
    else if(var == 9) then
      b(c) = b(c) + (P12 + VAR2_12 + VAR1w_12 + VAR2w_12)*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*1.5*(n1*n1+n2*n2) + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 
!==============================================================================================================================!
!==============================================================================================================================!
!------- uw stress
    else if(var == 10) then
      b(c) = b(c) + (P13 + VAR2_13 + VAR1w_13 + VAR2w_13)*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*1.5*(n1*n1+n3*n3) + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 

!==============================================================================================================================!
!==============================================================================================================================!
!------- vw stress
    else if(var == 11) then
      b(c) = b(c) + (P23 + VAR2_23 + VAR1w_23 + VAR2w_23)*volume(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*Eps%n(c)/Kin%n(c)+C1W*fw*Eps%n(c)/Kin%n(c)*1.5*(n2*n2+n3*n3) + &
                      fss*Eps%n(c)/Kin%n(c))*volume(c) 
!==============================================================================================================================!
!==============================================================================================================================!
!------- Eps eq.
    else if(var == 13) then
      Feps = 1.0 - ((Ce2-1.4)/Ce2)*exp(-(Ret/6.0)**2.0)
      Eps1 = 1.44*Pk(c)*Eps%n(c)/Kin%n(c)
      Eps2 = Ce2*Feps*Eps%n(c)/Kin%n(c)
      b(c) = b(c) + max(Eps1 + Diss1(c),0.0)*volume(c) 
     
      A % val(A % dia(c)) =  A % val(A % dia(c)) + Eps2*volume(c)
    end if
  end do

  do c = 1, NC
    VAR6x(c) = sqrt(0.5*(uu%n(c)+vv%n(c)+ww%n(c)))    
  end do 

  if(var == 13) then 
    call GraPhi(VAR6x,1,VAR7x,.TRUE.)             ! dK/dx
    call GraPhi(VAR6x,2,VAR7y,.TRUE.)             ! dK/dy
    call GraPhi(VAR6x,3,VAR7z,.TRUE.)             ! dK/dz
    do c = 1, NC
      Ret= (Kin % n(c)**2)/(VISc*Eps % n(c)+tiny)
      Feps = 1.0 - ((Ce2-1.4)/Ce2)*exp(-(Ret/6.0)**2.0)
      b(c) = b(c) + (Ce2*Feps*Eps%n(c)/Kin%n(c)*(VISc*(VAR7x(c)*VAR7x(c)+&
                    VAR7y(c)*VAR7y(c)+VAR7z(c)*VAR7z(c))))*volume(c)
    end do
  end if

  if(var == 13) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

!---- Calculate a values of dissipation  on wall
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          Eps%n(c2) = VISc*(VAR7x(c1)*VAR7x(c1)+VAR7y(c1)*VAR7y(c1)+VAR7z(c1)*VAR7z(c1)) 
        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if
  deallocate (Diss1)
  RETURN
  end subroutine SourcesHJ   


