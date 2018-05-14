!==============================================================================!
  subroutine Source_Ebm(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for 'EBM'.                                                     !
!   Following paper "Recent progress in the development of                     !
!   the Elliptic Blending Reynolds-stress model"                               !
!   DOI: doi.org/10.1016/j.ijheatfluidflow.2014.09.002                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Grad_Mod
  use Work_Mod, only: n1n1                  => r_cell_01, &
                      n2n2                  => r_cell_02, &
                      n3n3                  => r_cell_03, &
                      n1n2                  => r_cell_04, &
                      n1n3                  => r_cell_05, &
                      n2n3                  => r_cell_06, &
                      mag_f22               => r_cell_07, &
                      eps_2_k               => r_cell_08, &
                      alpha3                => r_cell_09, &
                      b11                   => r_cell_10, &
                      b22                   => r_cell_11, &
                      b33                   => r_cell_12, &
                      b12                   => r_cell_13, &
                      b13                   => r_cell_14, &
                      b23                   => r_cell_15, &
                      s11                   => r_cell_16, &
                      s22                   => r_cell_17, &
                      s33                   => r_cell_18, &
                      s12                   => r_cell_19, &
                      s13                   => r_cell_20, &
                      s23                   => r_cell_21, &
                      v12                   => r_cell_22, &
                      v13                   => r_cell_23, &
                      v23                   => r_cell_24, &
                      b_kl_b_kl_sq          => r_cell_25, &
                      b_lm_s_lm             => r_cell_26, &
                      u_k_u_l_n_k_n_l       => r_cell_27, &
                      term_c3_1             => r_cell_28
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2
  real    :: prod_and_coriolis, phi_hom, phi_wall, Esor
  real    :: c_1e1
  real    :: stress
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

  ! no need to compute this for EPS -> can be improved
  kin % n(1:) = max(0.5*(uu % n(1:) + vv % n(1:) + ww % n(1:)), TINY)
  ! P_k = 0.5 P_ii = - u_i u_k dU_i/dx_k
  p_kin(1:) = max(-( uu % n(1:) * u % x(1:)  &
                   + uv % n(1:) * u % y(1:)  &
                   + uw % n(1:) * u % z(1:)  &
                   + uv % n(1:) * v % x(1:)  &
                   + vv % n(1:) * v % y(1:)  &
                   + vw % n(1:) * v % z(1:)  &
                   + uw % n(1:) * w % x(1:)  &
                   + vw % n(1:) * w % y(1:)  &
                   + ww % n(1:) * w % z(1:)), TINY)

  ! |df22/x_j|
  mag_f22(1:) = max( f22 % x(1:)**2. + f22 % y(1:)**2. + f22 % z(1:)**2., TINY )

  ! formula C.9 (n_i never appears individually, only as n_i * n_j)
  n1n1(1:) = f22 % x(1:)**2.         / mag_f22(1:)
  n2n2(1:) = f22 % y(1:)**2.         / mag_f22(1:)
  n3n3(1:) = f22 % z(1:)**2.         / mag_f22(1:)
  n1n2(1:) = f22 % x(1:)*f22 % y(1:) / mag_f22(1:)
  n1n3(1:) = f22 % x(1:)*f22 % z(1:) / mag_f22(1:)
  n2n3(1:) = f22 % y(1:)*f22 % z(1:) / mag_f22(1:)

  ! frequently used expressions
  alpha3(1:)  = f22 % n(1:)**3.
  eps_2_k(1:) = eps % n(1:) / kin % n(1:)

  ! formula C.4
  b11(1:) = uu % n(1:)/(2.*kin % n(1:)) - ONE_THIRD 
  b22(1:) = vv % n(1:)/(2.*kin % n(1:)) - ONE_THIRD
  b33(1:) = ww % n(1:)/(2.*kin % n(1:)) - ONE_THIRD
  b12(1:) = uv % n(1:)/(2.*kin % n(1:))   
  b13(1:) = uw % n(1:)/(2.*kin % n(1:))    
  b23(1:) = vw % n(1:)/(2.*kin % n(1:))    

  ! formula C.5
  s11(1:) = u % x(1:) 
  s22(1:) = v % y(1:) 
  s33(1:) = w % z(1:) 
  s12(1:) = 0.5*(u % y(1:) + v % x(1:))
  s13(1:) = 0.5*(u % z(1:) + w % x(1:))
  s23(1:) = 0.5*(v % z(1:) + w % y(1:))

  ! formula C.6
  v12(1:) = 0.5*(u % y(1:) - v % x(1:)) - omega_z
  !v21 = -v12
  v13(1:) = 0.5*(u % z(1:) - w % x(1:)) + omega_y
  !v31 = -v13
  v23(1:) = 0.5*(v % z(1:) - w % y(1:)) - omega_x
  !v32 = -v23

  ! (C.3 1st term without "-" and *b_ij)
  term_c3_1(1:) = g1*eps % n(1:) + g1_star*p_kin(1:)

  ! for formula C.3 (b_kl_b_kl never appears without sqrt)
  b_kl_b_kl_sq(1:) = sqrt(b11(1:)**2. + b22(1:)**2. + b33(1:)**2. &
    + 2.*(b12(1:)**2. + b13(1:)**2. + b23(1:)**2.))

  ! for formula C.3
  b_lm_s_lm(1:) = b11(1:)*s11(1:) + b22(1:)*s22(1:) + b33(1:)*s33(1:) &
    + 2.*(b12(1:)*s12(1:) + b13(1:)*s13(1:) + b23(1:)*s23(1:))

  ! for formula C.7
  u_k_u_l_n_k_n_l(1:) = uu % n(1:)*n1n1(1:) + vv % n(1:)*n2n2(1:) &
    + ww % n(1:)*n3n3(1:) + 2.*uv % n(1:)*n1n2(1:) + 2.*uw % n(1:)*n1n3(1:) &
    + 2.*vw % n(1:)*n2n3(1:)

  do  c = 1, grid % n_cells
!------------------------------------------------------------------------------!
!   uu stress
!------------------------------------------------------------------------------!
    if (name_phi == 'UU') then
      ! useful vars
      stress = max(uu % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        2.*uu % n(c)*n1n1(c) + 2.*uv % n(c)*n1n2(c) + 2.*uw % n(c)*n1n3(c) &
        - 0.5*u_k_u_l_n_k_n_l(c)*(n1n1(c) + 1.) &
        - stress) ! this term is substracted in A later
      
      ! formula C.3 
      phi_hom = term_c3_1(c) * ONE_THIRD +  &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s11(c) + &
        (g4*(2.*(b11(c)*s11(c)+b12(c)*s12(c)+b13(c)*s13(c)) &
          - TWO_THIRDS*b_lm_s_lm(c)) + &
         g5*(2.*(              b12(c)*v12(c)+b13(c)*v13(c))))*kin % n(c)

      ! P_11 + G_11 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        -2.*(uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)) &
        -2.*omega_y*2.*uw % n(c) + 2.*omega_z*2.*uv % n(c) ! did not check

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3(c) * eps % n(c) / stress
!------------------------------------------------------------------------------!
!   vv stress
!------------------------------------------------------------------------------!
    elseif (name_phi == 'VV') then

      ! useful vars
      stress = max(vv % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        2.*uv % n(c)*n1n2(c) + 2.*vv % n(c)*n2n2(c) + 2.*vw % n(c)*n2n3(c)&
        - 0.5*u_k_u_l_n_k_n_l(c)*(n2n2(c) + 1.) &
        - stress) ! this term is substracted in A later

      ! formula C.3 
      phi_hom = term_c3_1(c) * ONE_THIRD +  &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s22(c) + &
        (g4*(2.*( b12(c)*s12(c)+b22(c)*s22(c)+b23(c)*s23(c)) &
          - TWO_THIRDS*b_lm_s_lm(c)) + &
         g5*(2.*(-b12(c)*v12(c)+              b23(c)*v23(c))))*kin % n(c)

      ! P_22 + G_22 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        -2.*(uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)) &
        -2.*omega_x*2.*vw % n(c) + 2.*omega_z*2.*uw % n(c)! did not check

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3(c) * eps % n(c) / stress
!------------------------------------------------------------------------------!
!   ww stress
!------------------------------------------------------------------------------!
    elseif (name_phi == 'WW') then

      ! useful vars
      stress = max(ww % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        2.*uw % n(c)*n1n3(c) + 2.*vw % n(c)*n2n3(c) + 2.*ww % n(c)*n3n3(c)&
        - 0.5*u_k_u_l_n_k_n_l(c)*(n3n3(c) + 1.) &
        - stress) ! this term is substracted in A later

      ! formula C.3 
      phi_hom = term_c3_1(c) * ONE_THIRD +  &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s33(c) + &
        (g4*(2.*( b13(c)*s13(c)+b23(c)*s23(c)+b33(c)*s33(c)) &
          - TWO_THIRDS*b_lm_s_lm(c)) + &
         g5*(2.*(-b13(c)*v13(c)-b23(c)*v23(c)        )))*kin % n(c)

      ! P_33 + G_33 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        -2.*(uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)) &
        -2.*omega_x*2.*vw % n(c) + 2.*omega_y*2.*uw % n(c)  ! did not check

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3(c) * eps % n(c) / stress
!------------------------------------------------------------------------------!
!   uv stress
!------------------------------------------------------------------------------!
    elseif (name_phi == 'UV') then

      ! useful vars
      stress = max(uv % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        uu % n(c)*n1n2(c) + uv % n(c)*n2n2(c) + uw % n(c)*n2n3(c) +  &
        uv % n(c)*n1n1(c) + vv % n(c)*n1n2(c) + vw % n(c)*n1n3(c)&
        - 0.5*u_k_u_l_n_k_n_l(c)*n1n2(c) &
        - stress) ! this term is substracted in A later

      ! formula C.3 
      phi_hom = &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s12(c) + &
        (g4*(b11(c)*s12(c)+b12(c)*s22(c)+b13(c)*s23(c)  + &
             b12(c)*s11(c)+b22(c)*s12(c)+b23(c)*s13(c)) + &
         g5*(b11(c)*v12(c)+              b13(c)*v23(c)  + &
                           b22(c)*v12(c)+b23(c)*v13(c)))*kin % n(c)

      ! P_12 + G_12 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        - uu % n(c)*v % x(c) - uv%n(c)*(v % y(c)+u % x(c)) - uw % n(c)*v % z(c)&
        - vv % n(c)*u % y(c) - vw % n(c)*u % z(c) &
        ! did not check
        + 2.*omega_x*uw%n(c)-2.*omega_y*vw%n(c)+2.*omega_z*(vv%n(c)-uu%n(c))
!------------------------------------------------------------------------------!
!   uw stress
!------------------------------------------------------------------------------!
    elseif (name_phi == 'UW') then

      ! useful vars
      stress = max(uw % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        uu % n(c)*n1n3(c) + uv % n(c)*n2n3(c) + uw % n(c)*n3n3(c)+ &
        uw % n(c)*n1n1(c)+ vw % n(c)*n1n2(c) + ww % n(c)*n1n3(c) &
        - 0.5*u_k_u_l_n_k_n_l(c)*n1n3(c) &
        - stress) ! this term is substracted in A later

      ! formula C.3 
      phi_hom = &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s13(c) + &
        (g4*(b11(c)*s13(c)+b12(c)*s23(c)+b13(c)*s33(c) + &
             b13(c)*s11(c)+b23(c)*s12(c)+b33(c)*s13(c)) + &
         g5*(b11(c)*v13(c)+b12(c)*v23(c)+       &
                           b23(c)*v12(c)+b33(c)*v13(c)))*kin % n(c)

      ! P_12 + G_12 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        ! did not check
        - 2.*omega_x*uv % n(c) - 2.*omega_y * &
        (ww % n(c) - uu % n(c)) + 2.*omega_z*vw % n(c)
!------------------------------------------------------------------------------!
!   vw stress
!------------------------------------------------------------------------------!
    elseif (name_phi == 'VW') then

      ! useful vars
      stress = max(vw % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k(c) * (  &
        vw % n(c) + uv % n(c)*n1n3(c) + vv % n(c)*n2n3(c) + vw % n(c)*n3n3(c) + &
        uw % n(c)*n1n2(c) + vw % n(c)*n2n2(c) + ww % n(c)*n2n3(c) &
        - 0.5*u_k_u_l_n_k_n_l(c)*n2n3(c) &
        - stress) ! this term is substracted in A later

      ! formula C.3 
      phi_hom = &
        (g3 - g3_star*b_kl_b_kl_sq(c))*kin % n(c)*s23(c) + &
        (g4*(b12(c)*s13(c)+b22(c)*s23(c)+b23(c)*s33(c) + &
             b13(c)*s12(c)+b23(c)*s22(c)+b33(c)*s23(c)) + &
         g5*(b12(c)*v13(c)+b22(c)*v23(c)+       &
             b13(c)*v12(c)+              b33(c)*v23(c)))*kin % n(c)

      ! P_12 + G_12 (formula C.1) + goes to right hand s. if > 0, on left if < 0
      prod_and_coriolis = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        ! did not check
        - 2.*omega_x*(vw % n(c)-ww % n(c)) &
        + 2.*omega_y*uv % n(c) - 2.*omega_z*uw % n(c)
!------------------------------------------------------------------------------!
!   repeating part for all stresses 
!------------------------------------------------------------------------------!
    ! formula C.1
    b(c) = b(c) + grid % vol(c) * ( & !
      max(prod_and_coriolis,0.) & ! P_ij + G_ij
      + (1. - alpha3(c))*phi_wall + alpha3(c)*phi_hom & ! C.2
      )
    ! left hand side
    A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * (  &
      + term_c3_1(c)/(2.*kin % n(c)) & ! from C.3 and C.4 1st terms
      - min(prod_and_coriolis,0.)/stress +  & ! (P_11 + G_11) / uu
      6.*(1. - alpha3(c))*eps_2_k(c) & ! C.11/uu
      )
!------------------------------------------------------------------------------!
!   eps 
!------------------------------------------------------------------------------!
    else if (name_phi == 'EPS') then
      Esor = grid % vol(c)/max(t_scale(c),TINY)
      c_1e1 = c_1e * (1. + 0.1*(1.-alpha3(c))*p_kin(c)/(eps % n(c)+TINY))
      b(c) = b(c) + c_1e1*density*p_kin(c)*Esor 

      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + c_2e*Esor*density
    end if 
  end do

  if (name_phi == 'EPS') then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Calculate values of dissipation on wall
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          eps % n(c2) = (viscosity/density)*&
            (uu % n(c1) + vv % n(c1) + ww % n(c1))/grid % wall_dist(c1)**2.
        end if ! end if of BC=wall
      end if   ! end if of c2<0
    end do
  end if ! name_phi == 'EPS'

  end subroutine