!==============================================================================!
  subroutine Update_Boundary_Values(grid)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: qx, qy, qz, Nx, Ny, Nz, Stot, CONeff, EBF, Yplus
  real    :: Prmol, beta, Uplus
!==============================================================================!

  Area  = 0.0
  Tflux = 0.0

  if(HOT==YES) then
    call GraPhi(grid, T % n,1, PHIx, .TRUE.)  ! dT/dx
    call GraPhi(grid, T % n,2, PHIy, .TRUE.)  ! dT/dy
    call GraPhi(grid, T % n,3, PHIz, .TRUE.)  ! dT/dz
  end if

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    !------------------------------------------------!
    !   Outflow (and inflow, if needed) boundaries   !
    !------------------------------------------------!

    ! On the boundary perform the extrapolation
    if(c2  < 0) then

      ! Extrapolate velocities on the outflow boundary 
      ! SYMMETRY is intentionally not treated here because I wanted to
      ! be sure that is handled only via graPHI and NewUVW functions)
      if( TypeBC(c2) == OUTFLOW.or.TypeBC(c2) == PRESSURE ) then
        U % n(c2) = U % n(c1)
        V % n(c2) = V % n(c1)
        W % n(c2) = W % n(c1)
        if(HOT==YES) T % n(c2) = T % n(c1)
      end if

      if( TypeBC(c2) == SYMMETRY ) then
        U % n(c2) = U % n(c1)
        V % n(c2) = V % n(c1)
        W % n(c2) = W % n(c1)
        if(HOT==YES) T % n(c2) = T % n(c1)
      end if

      ! Spalart Allmaras
      if(SIMULA == SPA_ALL .or. SIMULA == DES_SPA) then
        if ( TypeBC(c2) == OUTFLOW.or. TypeBC(c2) == CONVECT.or.&
        TypeBC(c2) == PRESSURE ) then
          VIS % n(c2) = VIS % n(c1) 
        end if
      end if

      if(SIMULA==EBM.or.SIMULA==HJ) then
        if ( TypeBC(c2) == WALL) then
          uu % n(c2) = 0.0 
          vv % n(c2) = 0.0 
          ww % n(c2) = 0.0 
          uv % n(c2) = 0.0 
          uw % n(c2) = 0.0 
          vw % n(c2) = 0.0 
          KIN % n(c2) = 0.0 
          if(SIMULA==EBM) f22% n(c2) = 0.0 
        end if
      end if

      ! k-epsilon-v^2
      if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
        if(TypeBC(c2) == OUTFLOW.or.TypeBC(c2) == CONVECT.or.&
        TypeBC(c2) == PRESSURE) then
          Kin % n(c2) = Kin % n(c1)
          Eps % n(c2) = Eps % n(c1)
          v_2 % n(c2) = v_2 % n(c1)
          f22 % n(c2) = f22 % n(c1)
        end if

        !  if (TypeBC(c2) == INFLOW) then
        !    f22 % n(c2) = f22 % n(c1)
        !  end if
      end if 

      ! k-epsilon
      if(SIMULA == K_EPS) then
        if(TypeBC(c2) == OUTFLOW.or.TypeBC(c2) == CONVECT.or.&
        TypeBC(c2) == PRESSURE.or.TypeBC(c2) == SYMMETRY) then
          Kin % n(c2) = Kin % n(c1)
          Eps % n(c2) = Eps % n(c1)
        end if
      end if 

      if(SIMULA==EBM.or.SIMULA==HJ) then
        if(TypeBC(c2) == OUTFLOW.or.TypeBC(c2) == CONVECT.or.&
           TypeBC(c2) == PRESSURE) then
          Kin % n(c2) = Kin % n(c1)
          Eps % n(c2) = Eps % n(c1)
          uu % n(c2) = uu % n(c1)
          vv % n(c2) = vv % n(c1)
          ww % n(c2) = ww % n(c1)
          uv % n(c2) = uv % n(c1)
          uw % n(c2) = uw % n(c1)
          vw % n(c2) = vw % n(c1)
          if(SIMULA==EBM) f22 % n(c2) = f22 % n(c1)
        end if
      end if 

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(HOT==YES) then
        Prt = 0.9
        if(SIMULA/=LES.or.SIMULA/=DNS) then
          Prt = 1.0                           &
              / ( 0.5882 + 0.228*(VISt(c1)    &
              / (VISc+1.0e-12)) - 0.0441      &
              * (VISt(c1)/(VISc+1.0e-12))**2  &
              * (1.0 - exp(-5.165*( VISc/(VISt(c1)+1.0e-12) ))))
        end if
        Stot = sqrt(  grid % sx(s)*grid % sx(s)  &
                    + grid % sy(s)*grid % sy(s)  &
                    + grid % sz(s)*grid % sz(s))
        Nx = grid % sx(s)/Stot
        Ny = grid % sy(s)/Stot
        Nz = grid % sz(s)/Stot
        qx = T % q(c2) * Nx 
        qy = T % q(c2) * Ny
        qz = T % q(c2) * Nz
        CONeff = CONc(material(c1))                 &
                + CAPc(material(c1))*VISt(c1)/Prt
        if(SIMULA==ZETA.or.SIMULA==K_EPS) then
          Yplus = max(Cmu25 * sqrt(Kin%n(c1)) * WallDs(c1)/VISc,0.12)
          Uplus = log(Yplus*Elog) / (kappa + TINY) + TINY
          Prmol = VISc / CONc(material(c1))
          beta = 9.24 * ((Prmol/Prt)**0.75 - 1.0)  &
                     * (1.0 + 0.28 * exp(-0.007*Prmol/Prt))
          EBF = 0.01 * (Prmol*Yplus)**4.0          &
                     / (1.0 + 5.0 * Prmol**3 * Yplus) + TINY
          CONwall(c1) = Yplus * VISc * CAPc(material(c1))   &
                      / (Yplus * Prmol * exp(-1.0 * EBF)    &
                      + (Uplus + beta) * Prt * exp(-1.0 / EBF) + TINY)
          if(TypeBC(c2) == WALLFL) then
            T% n(c2) = T % n(c1) + Prt / CAPc(material(c1))    &
                     * (  qx * grid % dx(s)   &
                        + qy * grid % dy(s)   &
                        + qz * grid % dz(s))  &
                     / CONwall(c1)
            Tflux = T % q(c2)
          else if(TypeBC(c2) == WALL) then
            T % q(c2) = ( T % n(c2) - T % n(c1) ) * CONeff     &
                      / (  Nx * grid % dx(s)  &
                         + Ny * grid % dy(s)  &
                         + Nz * grid % dz(s) )
            Tflux = T % q(c2)
          end if
        else
          if(TypeBC(c2) == WALLFL) then
            T% n(c2) = T % n(c1) + Prt / (CONeff + TINY)      &
                    * (  qx * grid % dx(s)  &
                       + qy * grid % dy(s)  &
                       + qz * grid % dz(s) ) 
            Tflux = T % q(c2) 
          else if(TypeBC(c2) == WALL) then
            T % q(c2) = ( T % n(c2) - T % n(c1) ) * CONeff     &
                      / (  Nx * grid % dx(s)  &
                         + Ny * grid % dy(s)  &
                         + Nz * grid % dz(s) )
            Tflux = T % q(c2) 
          end if
        end if
      end if

      ! Below part is coded for pipe flow with constant wall temperature    
      !
      !      if( (HOT==YES) .and. (TypeBC(c2) == WALL .or.  &
      !           TypeBC(c2) == WALLFL) ) then
      !        WallDs(c1) = (  grid % dx(s)**2  &
      !                      + grid % dy(s)**2  &
      !                      + grid % dz(s)**2)**0.5
      !        Stot = sqrt(  grid % sx(s)*grid % sx(s)  &
      !                    + grid % sy(s)*grid % sy(s)  &
      !                    + grid % sz(s)*grid % sz(s))
      !        if(Ynd(c1) < 3) then
      !          Tflux = Tflux + CONc(material(c1))*abs(PHIz(c1)) * Stot
      !        else 
      !          Tflux = Tflux + CONeff*(T % n(c1))  &
      !                / WallDs(c1) * Stot
      !        end if
      !        write(*,*) (T % n(c2) - T % n(c1))/WallDs(c1), PHIz(c1)
      !        Area  = Area  + Stot
      !      end if

      !---------------------!
      !   Copy boundaries   !
      !---------------------!
      if(CopyC(c2) /= 0) then
        U % n(c2) = U % n(CopyC(c2))
        V % n(c2) = V % n(CopyC(c2))
        W % n(c2) = W % n(CopyC(c2))
        if(HOT==YES)        T % n(c2) = T % n(CopyC(c2))
        if(SIMULA==SPA_ALL .or. SIMULA == DES_SPA) &
        VIS % n(c2) = VIS % n(CopyC(c2)) 
        if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
          Kin % n(c2) = Kin % n(CopyC(c2))
          Eps % n(c2) = Eps % n(CopyC(c2))
          v_2 % n(c2) = v_2 % n(CopyC(c2))
          f22 % n(c2) = f22 % n(CopyC(c2))
        end if ! K_EPS_VV
        if(SIMULA==K_EPS) then
          Kin % n(c2) = Kin % n(CopyC(c2))
          Eps % n(c2) = Eps % n(CopyC(c2))
        end if ! K_EPS
        if(SIMULA==EBM.or.SIMULA==HJ) then
          Kin % n(c2) = Kin % n(CopyC(c2))
          Eps % n(c2) = Eps % n(CopyC(c2))
          uu % n(c2) = uu % n(CopyC(c2))
          vv % n(c2) = vv % n(CopyC(c2))
          ww % n(c2) = ww % n(CopyC(c2))
          uv % n(c2) = uv % n(CopyC(c2))
          uw % n(c2) = uw % n(CopyC(c2))
          vw % n(c2) = vw % n(CopyC(c2))
          f22 % n(c2) = f22 % n(CopyC(c2))
        end if ! EBM 
      end if
    end if
  end do

  !  if(HOT==YES) then 
  !    call GloSum(Tflux)
  !    call GloSum(Area)
  !    call wait
  !   
  !    Tflux  = Tflux/(Area+TINY)
  !    Qflux  = Tflux * Area/(4.0*AreaZ(1))
  !  end if 
  !  write(*,*) Tflux
  !  This is done for pipe flow
  !  In order to be consistent we did not use ideal d*pi but Area/L where
  !  Area is total wall surface of pipe and L is lenght of pipe
  !  AreaZ(1) is surface in direction of to the flow stream

  end subroutine
