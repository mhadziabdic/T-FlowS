!==============================================================================!
  subroutine Source_Ebm(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for EBM.                                                       !  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use grid_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, i
  real    :: Prod, Diss, PHI_hom, PHI_wall, mag, PHI_tot, Esor
  real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32
  real    :: b11, b22, b33, b12, b13, b21, b31, b23, b32
  real    :: S11, S22, S33, S12, S13, S21, S31, S23, S32
  real    :: V11, V22, V33, V12, V13, V21, V31, V23, V32
  real    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uzy, Uzx, Uyx               
  real    :: n1, n2, n3, b_mn_b_mn, b_lk_s_lk, uiujn, Ce11, uu_nn 
  real    :: Diss_wall, Diss_hom, r13, r23
  real,allocatable :: Diss1(:)
!==============================================================================!

  call GraPhi(f22 % n, 1, VAR2x, .TRUE.)             ! df22/dx
  call GraPhi(f22 % n, 2, VAR2y, .TRUE.)             ! df22/dy
  call GraPhi(f22 % n, 3, VAR2z, .TRUE.)             ! df22/dz

  call Time_And_Length_Scale(grid)

  r13 = 1.0/3.0 
  do  c = 1, grid % n_cells
    Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
    Pk(c)      = max(-(uu % n(c)*Ux(c) + uv % n(c)*Uy(c) + uw % n(c)*Uz(c) +&
                   uv % n(c)*Vx(c) + vv % n(c)*Vy(c) + vw % n(c)*Vz(c) +&
                   uw % n(c)*Wx(c) + vw % n(c)*Wy(c) + ww % n(c)*Wz(c)),1.0e-12)                
  
    mag = max(1.0e-12,sqrt(VAR2x(c)*VAR2x(c)+VAR2y(c)*VAR2y(c)+VAR2z(c)*VAR2z(c)))       
    n1 = VAR2x(c)/mag 
    n2 = VAR2y(c)/mag 
    n3 = VAR2z(c)/mag 

    b11 = uu % n(c)/(2.0*Kin % n(c)) - r13 
    b22 = vv % n(c)/(2.0*Kin % n(c)) - r13
    b33 = ww % n(c)/(2.0*Kin % n(c)) - r13
    b12 = uv % n(c)/(2.0*Kin % n(c))   
    b21 = b12
    b13 = uw % n(c)/(2.0*Kin % n(c))    
    b31 = b13
    b23 = vw % n(c)/(2.0*Kin % n(c))    
    b32 = b23
    
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

    b_mn_b_mn = b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12+b13*b13+b23*b23)
    b_lk_s_lk = b11*S11 + b22*S22 + b33*S33 + 2.0*(b12*S12+b13*S13+b23*S23)
    uu_nn     = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3 &
               + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3 &
               + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)


    ! uu stress
    if(name_phi == 'UU') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-uu % n(c)+2.0*uu % n(c)*n1*n1 + 2.0*uv % n(c)*n1*n2 + 2.0*uw % n(c)*n1*n3 &
                 - 0.5*(n1*n1+1.0)*uu_nn)

      PHI_hom = -g1*Eps % n(c)*(-r13) - g1_star*Pk(c)*(-r13) +  &
                 (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S11 + &
                 g4*Kin%n(c)*(2.0*(b11*S11+b12*S12+b13*S13)-2.0/3.0 * b_lk_s_lk) +&
                 g5*Kin%n(c)*(2.0*(b11*V11+b12*V12+b13*V13))


      Prod = -2.0*(uu%n(c)*Ux(c) + uv % n(c)*Uy(c) + uw % n(c)*Uz(c))  &
             -2.0*omegaY*2.0*uw%n(c) + 2.0*omegaZ*2.0*uv%n(c)   


      Diss_wall = uu % n(c)/Kin % n(c) * Eps % n(c) 
      Diss_hom  = 2.0/3.0 * Eps % n(c)

      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-Prod,0.0)/max(uu%n(c),1.0e-12) + (1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c) +&         
                      f22%n(c)*f22%n(c)*(Diss_hom/max(uu%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))&
                      +g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)

    ! vv stress
    else if(name_phi == 'VV') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-vv % n(c)+2.0*uv % n(c)*n2*n1 + 2.0*vv % n(c)*n2*n2 + 2.0*vw % n(c)*n2*n3 &
                 - 0.5*(n2*n2+1.0)*uu_nn)

      PHI_hom =  -g1*Eps % n(c)*(-r13) - g1_star*Pk(c)*(-r13) +  &
                  (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S22 + &
                 g4*Kin%n(c)*(2.0*(b21*S21+b22*S22+b23*S23)-2.0/3.0 * b_lk_s_lk) +&
                 g5*Kin%n(c)*(2.0*(b21*V21+b22*V22+b23*V23))

      Prod = -2.0*(uv % n(c)*Vx(c) + vv%n(c)*Vy(c) + vw % n(c)*Vz(c))  &
             +2.0*omegaX*2.0*vw%n(c) - 2.0*omegaZ*2.0*uw%n(c)   


      PHI_tot = (1.0-f22 % n(c)*f22 % n(c))*PHI_wall &
               + f22 % n(c)*f22 % n(c)*PHI_hom 

      Diss_wall = vv % n(c)/Kin % n(c) * Eps % n(c) 
      Diss_hom  = 2.0/3.0 * Eps % n(c)

      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-Prod,0.0)/max(vv%n(c),1.0e-12) + (1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c) &         
                    + f22%n(c)*f22%n(c)*(Diss_hom/max(vv%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))&
                    +g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)

    ! ww stress
    else if(name_phi == 'WW') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-ww % n(c)+2.0*uw % n(c)*n3*n1 + 2.0*vw % n(c)*n3*n2 + 2.0*ww % n(c)*n3*n3 &
                 - 0.5*(n3*n3+1.0)*uu_nn)

      PHI_hom = -g1*Eps % n(c)*(-r13) - g1_star*Pk(c)*(-r13) +  &
                 (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S33 + &
                 g4*Kin%n(c)*(2.0*(b31*S31+b32*S32+b33*S33)-2.0/3.0 * b_lk_s_lk) +&
                 g5*Kin%n(c)*(2.0*(b31*V31+b32*V32+b33*V33))

      Prod = -2.0*(uw % n(c)*Wx(c) + vw % n(c)*Wy(c) + ww%n(c)*Wz(c))  &
             -2.0*omegaX*2.0*vw%n(c) + 2.0*omegaY*2.0*uw%n(c) 

      PHI_tot = (1.0-f22 % n(c)*f22 % n(c))*PHI_wall &
               + f22 % n(c)*f22 % n(c)*PHI_hom 

      Diss_wall = vv % n(c)/Kin % n(c) * Eps % n(c) 
      Diss_hom  = 2.0/3.0 * Eps % n(c)

      b(c) = b(c) + (max(Prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-Prod,0.0)/max(ww%n(c),1.0e-12)+(1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c)          &
                    + f22%n(c)*f22%n(c)*(Diss_hom/max(ww%n(c),1.0e-12)+g1*Eps%n(c)/(2.0*Kin%n(c))&
                    +g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)

    ! uv stress
    else if(name_phi == 'UV') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-uv % n(c)+uu % n(c)*n2*n1 + uv % n(c)*n2*n2 + uw % n(c)*n2*n3 + &
                  uv % n(c)*n1*n1 + vv % n(c)*n1*n2 + vw % n(c)*n1*n3   &
                 - 0.5*n1*n2*uu_nn)

      PHI_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S12 + &
                 g4*Kin%n(c)*(b11*S21+b12*S22+b13*S23 + &
                              b21*S11+b22*S12+b23*S13) +&
                 g5*Kin%n(c)*(b11*V21+b12*V22+b13*V23 + &
                              b21*V11+b22*V12+b23*V13)
 
      Prod = -(uu % n(c)*Vx(c) + uw % n(c)*Vz(c) + uv%n(c)*(Vy(c)+Ux(c)) +&
               vv % n(c)*Uy(c) + vw % n(c)*Uz(c)) &
             +2.0*omegaX*uw%n(c) - 2.0*omegaY*vw%n(c) + 2.0*omegaZ*(vv%n(c)-uu%n(c))  

      PHI_tot = (1.0-f22 % n(c)*f22 % n(c))*PHI_wall &
               + f22 % n(c)*f22 % n(c)*PHI_hom 
      
      Diss_wall = uv % n(c)/Kin % n(c) * Eps % n(c) 
 
      b(c) = b(c) + (Prod + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c)&
                    + f22%n(c)*f22%n(c)*( &
                    +g1*Eps%n(c)/(2.0*Kin%n(c))+g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)


    ! uw stress
    else if(name_phi == 'UW') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-uw % n(c)+uu % n(c)*n3*n1 + uv % n(c)*n3*n2 + uw % n(c)*n3*n3 + &
                  uw % n(c)*n1*n1 + vw % n(c)*n1*n2 + ww % n(c)*n1*n3   &
                 - 0.5*n1*n3*uu_nn)

      PHI_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S13 + &
                 g4*Kin%n(c)*(b11*S31+b12*S32+b13*S33 + &
                              b31*S11+b32*S12+b33*S13) +&
                 g5*Kin%n(c)*(b11*V31+b12*V32+b13*V33 + &
                              b31*V11+b32*V12+b33*V13)

      Prod = -(uu % n(c)*Wx(c) + uv % n(c)*Wy(c)+ uw%n(c)*(Wz(c)+Ux(c)) +&
               vw % n(c)*Uy(c) + ww % n(c)*Uz(c)) & 
             -2.0*omegaX*uv%n(c)-2.0*omegaY*(ww%n(c)-uu%n(c))+2.0*omegaZ*vw%n(c)   

      PHI_tot = (1.0-f22 % n(c)*f22 % n(c))*PHI_wall &
               + f22 % n(c)*f22 % n(c)*PHI_hom 

      Diss_wall = uw % n(c)/Kin % n(c) * Eps % n(c) 

      b(c) = b(c) + (Prod + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c)&           
                    + f22%n(c)*f22%n(c)*(&
                    +g1*Eps%n(c)/(2.0*Kin%n(c))+g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)

    ! vw stress
    else if(name_phi == 'VW') then
      PHI_wall = -5.0*Eps % n(c)/Kin % n(c) * &
                 (-vw % n(c)+uv % n(c)*n3*n1 + vv % n(c)*n3*n2 + vw % n(c)*n3*n3 + &
                             uw % n(c)*n2*n1 + vw % n(c)*n2*n2 + ww % n(c)*n2*n3   &
                 - 0.5*n2*n3*uu_nn)

      PHI_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*Kin % n(c)*S23 + &
                 g4*Kin%n(c)*(b21*S31+b22*S32+b23*S33 + &
                              b31*S21+b32*S22+b33*S23) +&
                 g5*Kin%n(c)*(b21*V31+b22*V32+b23*V33 + &
                              b31*V21+b32*V22+b33*V23)

      Prod = -(uv % n(c)*Wx(c) + vv % n(c)*Wy(c)+ vw%n(c)*(Wz(c)+Vy(c))+&
               uw % n(c)*Vx(c) + ww % n(c)*Vz(c))  &
             -2.0*omegaX*(vw%n(c)-ww%n(c))+2.0*omegaY*uv%n(c)-2.0*omegaZ*uw%n(c)   

      PHI_tot = (1.0-f22 % n(c)*f22 % n(c))*PHI_wall &
               + f22 % n(c)*f22 % n(c)*PHI_hom 

      Diss = (1.0 - f22 % n(c)*f22 % n(c)) * vw % n(c)/Kin % n(c) * Eps % n(c) 

      b(c) = b(c) + (Prod + (1.0-f22 % n(c)*f22 % n(c))*PHI_wall +&
             f22 % n(c)*f22 % n(c)*(PHI_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*Eps%n(c)/Kin%n(c) &           
                    + f22%n(c)*f22%n(c)*(&
                    +g1*Eps%n(c)/(2.0*Kin%n(c))+g1_star*Pk(c)/(2.0*Kin%n(c))))*grid % vol(c)

    ! Eps equation
    else if(name_phi == 'EPS') then
      Esor = grid % vol(c)/max(Tsc(c),1.0e-12)
!
!     Ce11 = Ce1*(1.0 + 0.03*(1. - f22%n(c)*f22%n(c))*sqrt(Kin%n(c)/max(uu_nn,1.e-12)))  
!
      Ce11 = Ce1*(1.0 + 0.1*(1.0-f22%n(c)**3.0)*Pk(c)/Eps%n(c))  
      b(c) = b(c) + Ce11*Pk(c)*Esor 

      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + Ce2*Esor*DENc(material(c))
    end if
  end do

  if(name_phi == 'EPS') then
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)

      ! Calculate values of dissipation on wall
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          Eps%n(c2) = VISc*(uu%n(c1)+vv%n(c1)+ww%n(c1))/WallDs(c1)**2
        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if

  end subroutine
