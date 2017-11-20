!==============================================================================!
  subroutine Compute_Sgs_Dynamic(grid)
!------------------------------------------------------------------------------!
!   Calculates Smagorinsky constant with dynamic procedure                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, j, cj
  real    :: u_a, v_a, w_a, uu_a, vv_a, ww_a, uv_a, uw_a, vw_a, vol_e
  real    :: m_11_a, m_22_a, m_33_a, m_12_a, m_13_a, m_23_a      
  real    :: l_11, l_22, l_33, l_12, l_13, l_23      
  real    :: m_11, m_22, m_33, m_12, m_13, m_23      
  real    :: m_dot_m, l_dot_m, Lg, Lf 
!==============================================================================!
!                                                                              !
!   C is derived from:    Lij_res = Lij_mod                                    !
!                                                                              !
!   res - resolved, mod - modeled                                              ! 
!                                                                              ! 
!   Lij_res = <UiUj> - <Ui><Uj>,                                               !
!                                                                              !
!   where <> denote test filter                                                !
!                                                                              !
!   Lij_mod = 2 * C * Mij, C = Csmag ** 2.0                                    !
!                                                                              !
!   where Mij is:  Mij = <delta**2.0>|<Sij>|<Sij> - <delta |Sij| Sij>          !
!                                                                              ! 
!   Finaly C is :                                                              ! 
!                                                                              !
!   C = 0.5 * Lij:Mij / Mij:Mij                                                !
!                                                                              !
!   Aij : Bij = A11 * B11 + A22 * B22 + A33 * B33                              !
!             + 2.0 * A12 * B12 + 2.0 A13 * B13 + 2.0 * A23 * B23              !   
!                                                                              !
!------------------------------------------------------------------------------!

  call Exchange(grid, U % n)
  call Exchange(grid, V % n)
  call Exchange(grid, W % n)
  call Exchange(grid, Shear)

  do c =1, grid % n_cells
    u_a   = 0.0
    v_a   = 0.0
    w_a   = 0.0
    vol_e = 0.0

    uu_a  = 0.0
    vv_a  = 0.0
    ww_a  = 0.0
    uv_a  = 0.0
    uw_a  = 0.0
    vw_a  = 0.0

    m_11_a = 0.0
    m_22_a = 0.0
    m_33_a = 0.0
    m_12_a = 0.0
    m_13_a = 0.0
    m_23_a = 0.0
  
    do j = A % row(c), A % row(c + 1) - 1
      cj = A % col(j) 
      if(cj /= c) then

        ! Test velocities
        u_a = u_a + grid % vol(cj) * U % n(cj)
        v_a = v_a + grid % vol(cj) * V % n(cj)
        w_a = w_a + grid % vol(cj) * W % n(cj)

        ! Test stresses
        uu_a = uu_a + grid % vol(cj) * U % n(cj) * U % n(cj)
        vv_a = vv_a + grid % vol(cj) * V % n(cj) * V % n(cj)
        ww_a = ww_a + grid % vol(cj) * W % n(cj) * W % n(cj)
        uv_a = uv_a + grid % vol(cj) * U % n(cj) * V % n(cj)
        uw_a = uw_a + grid % vol(cj) * U % n(cj) * W % n(cj)
        vw_a = vw_a + grid % vol(cj) * V % n(cj) * W % n(cj)

        ! Test Mija
        m_11_a = m_11_a + grid % vol(cj) * Shear(cj) * Ux(cj)
        m_22_a = m_22_a + grid % vol(cj) * Shear(cj) * Vy(cj)
        m_33_a = m_33_a + grid % vol(cj) * Shear(cj) * Wz(cj)
        m_12_a = m_12_a + grid % vol(cj) * Shear(cj)                &    
               * 0.5 * ( Uy(cj) + Vx(cj) ) 
        m_13_a = m_13_a + grid % vol(cj) * Shear(cj)                &
               * 0.5 * ( Uz(cj) + Wx(cj) )     
        m_23_a = m_23_a + grid % vol(cj) * Shear(cj)                &
               * 0.5 * ( Vz(cj) + Wy(cj) )

        ! Test volume 
        vol_e = vol_e + grid % vol(cj) 
      end if
    end do

    ! Take into account influence of central cell within test molecule
    vol_e = vol_e + grid % vol(c)

    u_a = u_a + grid % vol(c) * U % n(c)
    v_a = v_a + grid % vol(c) * V % n(c)
    w_a = w_a + grid % vol(c) * W % n(c)

    uu_a = uu_a + grid % vol(c) * U % n(c) * U % n(c)
    vv_a = vv_a + grid % vol(c) * V % n(c) * V % n(c)
    ww_a = ww_a + grid % vol(c) * W % n(c) * W % n(c)
    uv_a = uv_a + grid % vol(c) * U % n(c) * V % n(c)
    uw_a = uw_a + grid % vol(c) * U % n(c) * W % n(c)
    vw_a = vw_a + grid % vol(c) * V % n(c) * W % n(c)

    m_11_a = m_11_a + grid % vol(c) * Shear(c) * Ux(c)
    m_22_a = m_22_a + grid % vol(c) * Shear(c) * Vy(c)
    m_33_a = m_33_a + grid % vol(c) * Shear(c) * Wz(c)
    m_12_a = m_12_a + grid % vol(c) * Shear(c) * 0.5 * ( Uy(c) + Vx(c) ) 
    m_13_a = m_13_a + grid % vol(c) * Shear(c) * 0.5 * ( Uz(c) + Wx(c) )
    m_23_a = m_23_a + grid % vol(c) * Shear(c) * 0.5 * ( Vz(c) + Wy(c) )
    
    ! Now calculating test values
    U % filt(c) = u_a / vol_e
    V % filt(c) = v_a / vol_e
    W % filt(c) = w_a / vol_e

    UUf(c)  = uu_a / vol_e
    VVf(c)  = vv_a / vol_e
    WWf(c)  = ww_a / vol_e
    UVf(c)  = uv_a / vol_e
    UWf(c)  = uw_a / vol_e
    VWf(c)  = vw_a / vol_e
  
    M11f(c) = m_11_a / vol_e 
    M22f(c) = m_22_a / vol_e 
    M33f(c) = m_33_a / vol_e 
    M12f(c) = m_12_a / vol_e 
    M13f(c) = m_13_a / vol_e 
    M23f(c) = m_23_a / vol_e 
  end do

  call GraPhi(grid, U % filt, 1, Ux,.TRUE.)    ! dU/dx
  call GraPhi(grid, U % filt, 2, Uy,.TRUE.)    ! dU/dy
  call GraPhi(grid, U % filt, 3, Uz,.TRUE.)    ! dU/dz
  call GraPhi(grid, V % filt, 1, Vx,.TRUE.)    ! dV/dx
  call GraPhi(grid, V % filt, 2, Vy,.TRUE.)    ! dV/dy
  call GraPhi(grid, V % filt, 3, Vz,.TRUE.)    ! dV/dz
  call GraPhi(grid, W % filt, 1, Wx,.TRUE.)    ! dW/dx
  call GraPhi(grid, W % filt, 2, Wy,.TRUE.)    ! dW/dy
  call GraPhi(grid, W % filt, 3, Wz,.TRUE.)    ! dW/dz

  do c = 1, grid % n_cells
    Lg  = grid % vol(c)**ONE_THIRD
    Lf  = 2.0 * Lg

    ShearTest(c) = sqrt(2.0*(Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) +  &
                        0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) +           &
                        0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) +           &
                        0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c))))

    l_11 = UUf(c) - U % filt(c) * U % filt(c) 
    l_22 = VVf(c) - V % filt(c) * V % filt(c) 
    l_33 = WWf(c) - W % filt(c) * W % filt(c) 
    l_12 = UVf(c) - U % filt(c) * V % filt(c) 
    l_13 = UWf(c) - U % filt(c) * W % filt(c) 
    l_23 = VWf(c) - V % filt(c) * W % filt(c) 

    m_11 = Lf * Lf * ShearTest(c) * Ux(c)       &
         - Lg * Lg * M11f(c) 
    m_22 = Lf * Lf * ShearTest(c) * Vy(c)       &
         - Lg * Lg * M22f(c) 
    m_33 = Lf * Lf * ShearTest(c) * Wz(c)       &
         - Lg * Lg * M33f(c) 

    m_12 = Lf * Lf * ShearTest(c)                &
         * 0.5 * (Uy(c) + Vx(c)) - Lg * Lg * M12f(c)

    m_13 = Lf * Lf * ShearTest(c)                &
         * 0.5 * (Uz(c) + Wx(c)) - Lg * Lg * M13f(c) 

    m_23 = Lf * Lf * ShearTest(c)                &
         * 0.5 * (Vz(c) + Wy(c)) - Lg * Lg * M23f(c)  

    m_dot_m = m_11 * m_11 + m_22 * m_22 + m_33 * m_33   & 
            + 2.0 * m_12 * m_12                         &
            + 2.0 * m_13 * m_13                         &
            + 2.0 * m_23 * m_23 
 
    l_dot_m = l_11 * m_11 + l_22 * m_22 + l_33 * m_33   & 
            + 2.0 * l_12 * m_12                         &
            + 2.0 * l_13 * m_13                         &
            + 2.0 * l_23 * m_23

    Cdyn(c)  =  -0.5 * l_dot_m / (m_dot_m + TINY) 

    if(Cdyn(c) < 0.0) then
      Cdyn(c) = 0.0 
    else if(Cdyn(c) > 0.04) then
      Cdyn(c) = 0.04
    end if 
  end do

  end subroutine
