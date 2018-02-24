!==============================================================================!
  subroutine Source_Ebm(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for 'EBM'.                                                       !  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use allp_mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Grad_Mod
  use Work_Mod, only: f22_x => r_cell_01,  &
                      f22_y => r_cell_02,  &
                      f22_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, i
  real    :: prod, diss, phi_hom, phi_wall, mag, phi_tot, Esor
  real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32
  real    :: b11, b22, b33, b12, b13, b21, b31, b23, b32
  real    :: s11, s22, s33, s12, s13, s21, s31, s23, s32
  real    :: v11, v22, v33, v12, v13, v21, v31, v23, v32
  real    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uzy, Uzx, Uyx               
  real    :: n1, n2, n3, b_mn_b_mn, b_lk_s_lk, uiujn, Ce11, uu_nn 
  real    :: diss_wall, diss_hom, r13, r23
!==============================================================================!

  call Grad_Mod_For_Phi(grid, f22 % n, 1, f22_x, .true.)             ! df22/dx
  call Grad_Mod_For_Phi(grid, f22 % n, 2, f22_y, .true.)             ! df22/dy
  call Grad_Mod_For_Phi(grid, f22 % n, 3, f22_z, .true.)             ! df22/dz

  call Time_And_Length_Scale(grid)

  r13 = ONE_THIRD
  do  c = 1, grid % n_cells
    kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), 1.0e-12)
    p_kin(c)   = max(-(  uu % n(c) * u % x(c)  &
                       + uv % n(c) * u % y(c)  &
                       + uw % n(c) * u % z(c)  &
                       + uv % n(c) * v % x(c)  &
                       + vv % n(c) * v % y(c)  &
                       + vw % n(c) * v % z(c)  &
                       + uw % n(c) * w % x(c)  &
                       + vw % n(c) * w % y(c)  &
                       + ww % n(c) * w % z(c)), 1.0e-12)                
  
    mag = max(1.0e-12, sqrt(  f22_x(c)**2 + f22_y(c)**2 + f22_z(c)**2))       
    n1 = f22_x(c) / mag 
    n2 = f22_y(c) / mag 
    n3 = f22_z(c) / mag 

    b11 = uu % n(c)/(2.0*kin % n(c)) - r13 
    b22 = vv % n(c)/(2.0*kin % n(c)) - r13
    b33 = ww % n(c)/(2.0*kin % n(c)) - r13
    b12 = uv % n(c)/(2.0*kin % n(c))   
    b21 = b12
    b13 = uw % n(c)/(2.0*kin % n(c))    
    b31 = b13
    b23 = vw % n(c)/(2.0*kin % n(c))    
    b32 = b23
    
    s11 = u % x(c) 
    s22 = v % y(c) 
    s33 = w % z(c) 
    s12 = 0.5*(u % y(c)+v % x(c)) 
    s21 = s12
    s13 = 0.5*(u % z(c)+w % x(c)) 
    s31 = s13 
    s23 = 0.5*(v % z(c)+w % y(c)) 
    s32 = s23 

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 = 0.5*(u % y(c)-v % x(c)) - omega_z
    v21 = -v12 + omega_z
    v13 = 0.5*(u % z(c)-w % x(c)) + omega_y
    v31 = -v13 - omega_y
    v23 = 0.5*(v % z(c)-w % y(c)) - omega_x
    v32 = -v23 + omega_x

    b_mn_b_mn = b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12+b13*b13+b23*b23)
    b_lk_s_lk = b11*s11 + b22*s22 + b33*s33 + 2.0*(b12*s12+b13*s13+b23*s23)
    uu_nn     = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3 &
               + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3 &
               + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)


    ! uu stress
    if(name_phi == 'UU') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-uu % n(c)+2.0*uu % n(c)*n1*n1 + 2.0*uv % n(c)*n1*n2 + 2.0*uw % n(c)*n1*n3 &
                 - 0.5*(n1*n1+1.0)*uu_nn)

      phi_hom = -g1*eps % n(c)*(-r13) - g1_star*p_kin(c)*(-r13) +  &
                 (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s11 + &
                 g4*kin%n(c)*(2.0*(b11*s11+b12*s12+b13*s13)-2.0/3.0 * b_lk_s_lk) +&
                 g5*kin%n(c)*(2.0*(b11*v11+b12*v12+b13*v13))


      prod = -2.0*(uu%n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c))  &
             -2.0*omega_y*2.0*uw%n(c) + 2.0*omega_z*2.0*uv%n(c)   


      diss_wall = uu % n(c)/kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + (max(prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-prod,0.0)/max(uu%n(c),1.0e-12) + (1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c) +&         
                      f22%n(c)*f22%n(c)*(diss_hom/max(uu%n(c),1.0e-12)+g1*eps%n(c)/(2.0*kin%n(c))&
                      +g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)

    ! vv stress
    else if(name_phi == 'VV') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-vv % n(c)+2.0*uv % n(c)*n2*n1 + 2.0*vv % n(c)*n2*n2 + 2.0*vw % n(c)*n2*n3 &
                 - 0.5*(n2*n2+1.0)*uu_nn)

      phi_hom =  -g1*eps % n(c)*(-r13) - g1_star*p_kin(c)*(-r13) +  &
                  (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s22 + &
                 g4*kin%n(c)*(2.0*(b21*s21+b22*s22+b23*s23)-2.0/3.0 * b_lk_s_lk) +&
                 g5*kin%n(c)*(2.0*(b21*v21+b22*v22+b23*v23))

      prod = -2.0*(uv % n(c)*v % x(c) + vv%n(c)*v % y(c) + vw % n(c)*v % z(c))  &
             +2.0*omega_x*2.0*vw%n(c) - 2.0*omega_z*2.0*uw%n(c)   


      phi_tot = (1.0-f22 % n(c)*f22 % n(c))*phi_wall &
               + f22 % n(c)*f22 % n(c)*phi_hom 

      diss_wall = vv % n(c)/kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + (max(prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-prod,0.0)/max(vv%n(c),1.0e-12) + (1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c) &         
                    + f22%n(c)*f22%n(c)*(diss_hom/max(vv%n(c),1.0e-12)+g1*eps%n(c)/(2.0*kin%n(c))&
                    +g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)

    ! ww stress
    else if(name_phi == 'WW') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-ww % n(c)+2.0*uw % n(c)*n3*n1 + 2.0*vw % n(c)*n3*n2 + 2.0*ww % n(c)*n3*n3 &
                 - 0.5*(n3*n3+1.0)*uu_nn)

      phi_hom = -g1*eps % n(c)*(-r13) - g1_star*p_kin(c)*(-r13) +  &
                 (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s33 + &
                 g4*kin%n(c)*(2.0*(b31*s31+b32*s32+b33*s33)-2.0/3.0 * b_lk_s_lk) +&
                 g5*kin%n(c)*(2.0*(b31*v31+b32*v32+b33*v33))

      prod = -2.0*(uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww%n(c)*w % z(c))  &
             -2.0*omega_x*2.0*vw%n(c) + 2.0*omega_y*2.0*uw%n(c) 

      phi_tot = (1.0-f22 % n(c)*f22 % n(c))*phi_wall &
               + f22 % n(c)*f22 % n(c)*phi_hom 

      diss_wall = vv % n(c)/kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + (max(prod,0.0) + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + (max(-prod,0.0)/max(ww%n(c),1.0e-12)+(1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c)          &
                    + f22%n(c)*f22%n(c)*(diss_hom/max(ww%n(c),1.0e-12)+g1*eps%n(c)/(2.0*kin%n(c))&
                    +g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)

    ! uv stress
    else if(name_phi == 'UV') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-uv % n(c)+uu % n(c)*n2*n1 + uv % n(c)*n2*n2 + uw % n(c)*n2*n3 + &
                  uv % n(c)*n1*n1 + vv % n(c)*n1*n2 + vw % n(c)*n1*n3   &
                 - 0.5*n1*n2*uu_nn)

      phi_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s12 + &
                 g4*kin%n(c)*(b11*s21+b12*s22+b13*s23 + &
                              b21*s11+b22*s12+b23*s13) +&
                 g5*kin%n(c)*(b11*v21+b12*v22+b13*v23 + &
                              b21*v11+b22*v12+b23*v13)
 
      prod = -(uu % n(c)*v % x(c) + uw % n(c)*v % z(c) + uv%n(c)*(v % y(c)+u % x(c)) +&
               vv % n(c)*u % y(c) + vw % n(c)*u % z(c)) &
             +2.0*omega_x*uw%n(c) - 2.0*omega_y*vw%n(c) + 2.0*omega_z*(vv%n(c)-uu%n(c))  

      phi_tot = (1.0-f22 % n(c)*f22 % n(c))*phi_wall &
               + f22 % n(c)*f22 % n(c)*phi_hom 
      
      diss_wall = uv % n(c)/kin % n(c) * eps % n(c) 
 
      b(c) = b(c) + (prod + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c)&
                    + f22%n(c)*f22%n(c)*( &
                    +g1*eps%n(c)/(2.0*kin%n(c))+g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)


    ! uw stress
    else if(name_phi == 'UW') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-uw % n(c)+uu % n(c)*n3*n1 + uv % n(c)*n3*n2 + uw % n(c)*n3*n3 + &
                  uw % n(c)*n1*n1 + vw % n(c)*n1*n2 + ww % n(c)*n1*n3   &
                 - 0.5*n1*n3*uu_nn)

      phi_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s13 + &
                 g4*kin%n(c)*(b11*s31+b12*s32+b13*s33 + &
                              b31*s11+b32*s12+b33*s13) +&
                 g5*kin%n(c)*(b11*v31+b12*v32+b13*v33 + &
                              b31*v11+b32*v12+b33*v13)

      prod = -(uu % n(c)*w % x(c) + uv % n(c)*w % y(c)+ uw%n(c)*(w % z(c)+u % x(c)) +&
               vw % n(c)*u % y(c) + ww % n(c)*u % z(c)) & 
             -2.0*omega_x*uv%n(c)-2.0*omega_y*(ww%n(c)-uu%n(c))+2.0*omega_z*vw%n(c)   

      phi_tot = (1.0-f22 % n(c)*f22 % n(c))*phi_wall &
               + f22 % n(c)*f22 % n(c)*phi_hom 

      diss_wall = uw % n(c)/kin % n(c) * eps % n(c) 

      b(c) = b(c) + (prod + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c)&           
                    + f22%n(c)*f22%n(c)*(&
                    +g1*eps%n(c)/(2.0*kin%n(c))+g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)

    ! vw stress
    else if(name_phi == 'VW') then
      phi_wall = -5.0*eps % n(c)/kin % n(c) * &
                 (-vw % n(c)+uv % n(c)*n3*n1 + vv % n(c)*n3*n2 + vw % n(c)*n3*n3 + &
                             uw % n(c)*n2*n1 + vw % n(c)*n2*n2 + ww % n(c)*n2*n3   &
                 - 0.5*n2*n3*uu_nn)

      phi_hom = &
                 (g3-g3_star*sqrt(b_mn_b_mn))*kin % n(c)*s23 + &
                 g4*kin%n(c)*(b21*s31+b22*s32+b23*s33 + &
                              b31*s21+b32*s22+b33*s23) +&
                 g5*kin%n(c)*(b21*v31+b22*v32+b23*v33 + &
                              b31*v21+b32*v22+b33*v23)

      prod = -(uv % n(c)*w % x(c) + vv % n(c)*w % y(c)+ vw%n(c)*(w % z(c)+v % y(c))+&
               uw % n(c)*v % x(c) + ww % n(c)*v % z(c))  &
             -2.0*omega_x*(vw%n(c)-ww%n(c))+2.0*omega_y*uv%n(c)-2.0*omega_z*uw%n(c)   

      phi_tot = (1.0-f22 % n(c)*f22 % n(c))*phi_wall &
               + f22 % n(c)*f22 % n(c)*phi_hom 

      diss = (1.0 - f22 % n(c)*f22 % n(c)) * vw % n(c)/kin % n(c) * eps % n(c) 

      b(c) = b(c) + (prod + (1.0-f22 % n(c)*f22 % n(c))*phi_wall +&
             f22 % n(c)*f22 % n(c)*(phi_hom))*grid % vol(c)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + ((1.0-f22%n(c)*f22%n(c))*&
                      6.0*eps%n(c)/kin%n(c) &           
                    + f22%n(c)*f22%n(c)*(&
                    +g1*eps%n(c)/(2.0*kin%n(c))+g1_star*p_kin(c)/(2.0*kin%n(c))))*grid % vol(c)

    ! eps equation
    else if(name_phi == 'EPS') then
      Esor = grid % vol(c)/max(Tsc(c),1.0e-12)
!
!     Ce11 = Ce1*(1.0 + 0.03*(1. - f22%n(c)*f22%n(c))*sqrt(kin%n(c)/max(uu_nn,1.e-12)))  
!
      Ce11 = Ce1*(1.0 + 0.1*(1.0-f22%n(c)**3.0)*p_kin(c)/eps%n(c))  
      b(c) = b(c) + Ce11*p_kin(c)*Esor 

      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + Ce2*Esor*density
    end if
  end do

  if(name_phi == 'EPS') then
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)

      ! Calculate values of dissipation on wall
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then
          eps%n(c2) = viscosity*(uu%n(c1)+vv%n(c1)+ww%n(c1))/grid % wall_dist(c1)**2
        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if

  end subroutine
