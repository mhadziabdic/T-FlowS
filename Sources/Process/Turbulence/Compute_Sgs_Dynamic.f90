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
  integer           :: c, j, cj, m, n 
  real              :: Ua, Va, Wa, UUa, VVa, WWa, UVa, UWa, VWa, DE
  real              :: M11a, M22a, M33a, M12a, M13a, M23a      
  real              :: L11, L22, L33, L12, L13, L23      
  real              :: M11, M22, M33, M12, M13, M23      
  real              :: MdotM, LdotM, Lg, Lf 
  real              :: MinC, MaxC, Cplus_avr, Cminus_avr
  real              :: fun

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
    Ua   = 0.0
    Va   = 0.0
    Wa   = 0.0
    DE   = 0.0

    UUa  = 0.0
    VVa  = 0.0
    WWa  = 0.0
    UVa  = 0.0
    UWa  = 0.0
    VWa  = 0.0

    M11a = 0.0
    M22a = 0.0
    M33a = 0.0
    M12a = 0.0
    M13a = 0.0
    M23a = 0.0
  
    do j = A % row(c), A % row(c + 1) - 1
      cj = A % col(j) 
      if(cj /= c) then

        ! Test velocitys
        Ua = Ua + grid % vol(cj) * U % n(cj)
        Va = Va + grid % vol(cj) * V % n(cj)
        Wa = Wa + grid % vol(cj) * W % n(cj)

        ! Test stresses
        UUa = UUa + grid % vol(cj) * U % n(cj) * U % n(cj)
        VVa = VVa + grid % vol(cj) * V % n(cj) * V % n(cj)
        WWa = WWa + grid % vol(cj) * W % n(cj) * W % n(cj)
        UVa = UVa + grid % vol(cj) * U % n(cj) * V % n(cj)
        UWa = UWa + grid % vol(cj) * U % n(cj) * W % n(cj)
        VWa = VWa + grid % vol(cj) * V % n(cj) * W % n(cj)

        ! Test Mija
        M11a = M11a + grid % vol(cj) * Shear(cj) * Ux(cj)
        M22a = M22a + grid % vol(cj) * Shear(cj) * Vy(cj)
        M33a = M33a + grid % vol(cj) * Shear(cj) * Wz(cj)
        M12a = M12a + grid % vol(cj) * Shear(cj)                &    
               * 0.5 * ( Uy(cj) + Vx(cj) ) 
        M13a = M13a + grid % vol(cj) * Shear(cj)                &
               * 0.5 * ( Uz(cj) + Wx(cj) )     
        M23a = M23a + grid % vol(cj) * Shear(cj)                &
               * 0.5 * ( Vz(cj) + Wy(cj) )

        ! Test volume 
        DE = DE + grid % vol(cj) 
      end if
    end do

    ! Take into account influence of central cell within test molecule
    DE = DE + grid % vol(c)

    Ua = Ua + grid % vol(c) * U % n(c)
    Va = Va + grid % vol(c) * V % n(c)
    Wa = Wa + grid % vol(c) * W % n(c)

    UUa = UUa + grid % vol(c) * U % n(c) * U % n(c)
    VVa = VVa + grid % vol(c) * V % n(c) * V % n(c)
    WWa = WWa + grid % vol(c) * W % n(c) * W % n(c)
    UVa = UVa + grid % vol(c) * U % n(c) * V % n(c)
    UWa = UWa + grid % vol(c) * U % n(c) * W % n(c)
    VWa = VWa + grid % vol(c) * V % n(c) * W % n(c)

    M11a = M11a + grid % vol(c) * Shear(c) * Ux(c)
    M22a = M22a + grid % vol(c) * Shear(c) * Vy(c)
    M33a = M33a + grid % vol(c) * Shear(c) * Wz(c)
    M12a = M12a + grid % vol(c) * Shear(c) * 0.5 * ( Uy(c) + Vx(c) ) 
    M13a = M13a + grid % vol(c) * Shear(c) * 0.5 * ( Uz(c) + Wx(c) )
    M23a = M23a + grid % vol(c) * Shear(c) * 0.5 * ( Vz(c) + Wy(c) )
    
    ! Now calculating test values
    U % filt(c) = Ua / DE
    V % filt(c) = Va / DE
    W % filt(c) = Wa / DE

    UUf(c)      = UUa / DE
    VVf(c)      = VVa / DE
    WWf(c)      = WWa / DE
    UVf(c)      = UVa / DE
    UWf(c)      = UWa / DE
    VWf(c)      = VWa / DE
  
    M11f(c)     = M11a / DE 
    M22f(c)     = M22a / DE 
    M33f(c)     = M33a / DE 
    M12f(c)     = M12a / DE 
    M13f(c)     = M13a / DE 
    M23f(c)     = M23a / DE 
  end do

  call GraPhi(U % filt, 1, Ux,.TRUE.)    ! dU/dx
  call GraPhi(U % filt, 2, Uy,.TRUE.)    ! dU/dy
  call GraPhi(U % filt, 3, Uz,.TRUE.)    ! dU/dz
  call GraPhi(V % filt, 1, Vx,.TRUE.)    ! dV/dx
  call GraPhi(V % filt, 2, Vy,.TRUE.)    ! dV/dy
  call GraPhi(V % filt, 3, Vz,.TRUE.)    ! dV/dz
  call GraPhi(W % filt, 1, Wx,.TRUE.)    ! dW/dx
  call GraPhi(W % filt, 2, Wy,.TRUE.)    ! dW/dy
  call GraPhi(W % filt, 3, Wz,.TRUE.)    ! dW/dz

  do c = 1, grid % n_cells
    Lg  = grid % vol(c)**ONE_THIRD
    Lf  = 2.0 * Lg

    ShearTest(c) = sqrt(2.0*(Ux(c)*Ux(c) + Vy(c)*Vy(c) + Wz(c)*Wz(c) + &
             0.5*(Vz(c) + Wy(c))*(Vz(c) + Wy(c)) + &
             0.5*(Uz(c) + Wx(c))*(Uz(c) + Wx(c)) + &
             0.5*(Vx(c) + Uy(c))*(Vx(c) + Uy(c))))

    L11 = UUf(c) - U % filt(c) * U % filt(c) 
    L22 = VVf(c) - V % filt(c) * V % filt(c) 
    L33 = WWf(c) - W % filt(c) * W % filt(c) 
    L12 = UVf(c) - U % filt(c) * V % filt(c) 
    L13 = UWf(c) - U % filt(c) * W % filt(c) 
    L23 = VWf(c) - V % filt(c) * W % filt(c) 

    M11 = Lf * Lf * ShearTest(c) * Ux(c)       &
        - Lg * Lg * M11f(c) 
    M22 = Lf * Lf * ShearTest(c) * Vy(c)       &
        - Lg * Lg * M22f(c) 
    M33 = Lf * Lf * ShearTest(c) * Wz(c)       &
        - Lg * Lg * M33f(c) 

    M12 = Lf * Lf * ShearTest(c)                &
          * 0.5 * (Uy(c) + Vx(c)) - Lg * Lg * M12f(c)

    M13 = Lf * Lf * ShearTest(c)                &
          * 0.5 * (Uz(c) + Wx(c)) - Lg * Lg * M13f(c) 

    M23 = Lf * Lf * ShearTest(c)                &
          * 0.5 * (Vz(c) + Wy(c)) - Lg * Lg * M23f(c)  

    MdotM = M11 * M11 + M22 * M22 + M33 * M33    & 
          + 2.0 * M12 * M12                      &
          + 2.0 * M13 * M13                      &
          + 2.0 * M23 * M23 

    LdotM = L11 * M11 + L22 * M22 + L33 * M33    & 
          + 2.0 * L12 * M12                      &
          + 2.0 * L13 * M13                      &
          + 2.0 * L23 * M23

    Cdyn(c)  =  -0.5 * LdotM / (MdotM + TINY) 

    if(Cdyn(c) < 0.0) then
      Cdyn(c) = 0.0 
    else if(Cdyn(c) > 0.04) then
      Cdyn(c) = 0.04
    end if 
  end do

  end subroutine
