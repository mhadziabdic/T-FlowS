!==============================================================================!
  subroutine Compute_Sgs(grid)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for LES.                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Constants_Pro_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03           
!------------------------------------------------------------------------------!
!   Near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the near(c) is zero.                                   !
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig       !
!   of a simulation.                                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2 
  real    :: Nx, Ny, Nz
  real    :: Cs
  real    :: Stot, lf, UtauL, Uff 
  real    :: Utot, Unor, Utan, Apow, Bpow, nu, dely, yPlus 
  real    :: Nc2
!==============================================================================!
  
  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(BUOY==YES) then
    call GraPhi(grid, t % n, 1, t_x, .true.)  ! dT/dx
    call GraPhi(grid, t % n, 2, t_y, .true.)  ! dT/dy
    call GraPhi(grid, t % n, 3, t_z, .true.)  ! dT/dz
  end if 
 
  if(MODE == SMAG) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD

      ! if(near(c) /= 0) is needed for parallel version
      ! since the subdomains which do not "touch" wall
      ! has near(c) = 0. 
      if(near(c) /= 0) then
        Uff = sqrt( VISc  &
                    * sqrt(  U % n(near(c)) * U % n(near(c))    & 
                           + V % n(near(c)) * V % n(near(c))    &
                           + W % n(near(c)) * W % n(near(c)))   &
                   / (grid % wall_dist(near(c))+TINY) )
        yPlus = grid % wall_dist(c) * Uff / VISc
        Cs = Cs0 * (1.0-exp(-yPlus/25.0))
      else  
        Cs = Cs0
      end if
      vis_t(c) = DENc(material(c)) &
              * (lf*lf)           &          ! delta^2 
              * (Cs*Cs)           &          ! Cs^2   
              * Shear(c) 
    end do
  else if(MODE == DYN) then
    if(BUOY==YES) then  
      do c = 1, grid % n_cells
        lf = grid % vol(c)**ONE_THIRD  
        vis_t(c) = DENc(material(c))  &
                * (lf*lf)            &          ! delta^2 
                * Cdyn(c)            &          ! Cdynamic   
                * Shear(c) 
      end do
    else
      do c = 1, grid % n_cells 
        vis_t(c) = DENc(material(c))       &
                * (lf*lf)                 &          ! delta^2 
                * Cdyn(c)                 &          ! Cdynamic   
                * sqrt(Shear(c)*Shear(c)  &
                + 2.5*(grav_x*t_x(c) + grav_y*t_y(c) + grav_z*t_z(c)))  
      end do
    end if     
  else if(MODE == WALE) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD    
      vis_t(c) = DENc(material(c))  &
              * (lf*lf)            &          ! delta^2 
              * (0.5*0.5)          &          ! Cs^2   
              * WALEv(c) 
    end do
  end if

  if(BUOY==YES) then
    do c = 1, grid % n_cells
      Nc2 = -(  grav_x * t_x(c)   &
              + grav_y * t_y(c)   &
              + grav_z * t_z(c))  &
          / Tref
      Nc2 = max(0.0, Nc2) 
      vis_t(c) = vis_t(c) * sqrt(1.0 - min(2.5*Nc2/(Shear(c)**2), 1.0))
    end do
  end if

  !-------------------!
  !   Wall function   !
  !-------------------+--------------!
  !   Law of the wall:               !
  !                                  !
  !   u+ = yz+  for z+ < 11.81       !
  !                                  !
  !   and                            !
  !                                  !
  !   u+ = A(y+)^B   for y+ > 11.81  !
  !                                  !
  !   with: A = 8.3 and B = 1/7      !
  !                                  !
  !----------------------------------+----------!
  !   The procedure below should be activated   !
  !   only if wall function approach is used.   !
  !----------------.----------------------------! 
  do s = 1, grid % n_faces
    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s)

    if(c2  < 0) then 

      Stot = sqrt(grid % sx(s)*grid % sx(s) +  &
                  grid % sy(s)*grid % sy(s) +  &
                  grid % sz(s)*grid % sz(s))

      nx = grid % sx(s) / Stot 
      ny = grid % sy(s) / Stot 
      nz = grid % sz(s) / Stot 

      if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

        Utot = sqrt(  U % n(c1) * U % n(c1)     &
                    + V % n(c1) * V % n(c1)     & 
                    + W % n(c1) * W % n(c1)  )

        Unor = (   U % n(c1) * Nx  &
                 + V % n(c1) * Ny  &
                 + W % n(c1) * Nz )   

        if( abs(Utot) > abs(Unor) ) then
          Utan = sqrt(Utot * Utot - Unor * Unor)
        else
          Utan = TINY 
        end if

        Apow = 8.3
        Bpow = 1.0/7.0
        nu = VISc/DENc(material(c1))
        dely = grid % wall_dist(c1)

        ! Calculate UtauL
        UtauL = ( Utan/Apow * (nu/dely)**Bpow )                     &
          ** (1.0/(1.0+Bpow))

        ! Calculate TauWall 
        TauWall(c1) = VISc * Utan / dely 
 
        ! Calculate y+
        yPlus  = dely*UtauL/nu
        if(yPlus  >=  11.81) then
          ! This one is effective viscosity
          VISwall(c1) = DENc(material(c1))*UtauL*UtauL*dely/abs(Utan) 
        else 
          VISwall(c1) = VISc + fF(s)*vis_t(c1)+(1.0-fF(s))*vis_t(c2)
        endif
      end if  ! TypeBC(c2)==WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Exchange(grid, vis_t)
  call Exchange(grid, VISwall)

  end subroutine
