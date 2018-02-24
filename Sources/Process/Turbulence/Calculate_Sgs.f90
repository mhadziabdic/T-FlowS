!==============================================================================!
  subroutine Calculate_Sgs(grid)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03           
!------------------------------------------------------------------------------!
!   Near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                                   !
!   nearest_wall_cell(c) is calculated in NearWallCells.f90, only ones in the beginig       !
!   of a simulation.                                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2 
  real    :: nx, ny, nz
  real    :: Cs
  real    :: s_tot, lf, UtauL, Uff 
  real    :: u_tot, u_nor, u_tan, a_pow, b_pow, nu, dely, y_plus 
  real    :: Nc2
!==============================================================================!
  
  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(buoyancy == YES) then
    call GraPhi(grid, t % n, 1, t_x, .true.)  ! dT/dx
    call GraPhi(grid, t % n, 2, t_y, .true.)  ! dT/dy
    call GraPhi(grid, t % n, 3, t_z, .true.)  ! dT/dz
  end if 
 
  if(turbulence_model_variant == SMAGORINSKY) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD

      ! if(nearest_wall_cell(c) /= 0) is needed for parallel version
      ! since the subdomains which do not "touch" wall
      ! has nearest_wall_cell(c) = 0. 
      if(nearest_wall_cell(c) /= 0) then
        Uff = sqrt( viscosity  &
                    * sqrt(  u % n(nearest_wall_cell(c)) ** 2   & 
                           + v % n(nearest_wall_cell(c)) ** 2   &
                           + w % n(nearest_wall_cell(c)) ** 2)  &
                   / (grid % wall_dist(nearest_wall_cell(c))+TINY) )
        y_plus = grid % wall_dist(c) * Uff / viscosity
        Cs = Cs0 * (1.0-exp(-y_plus/25.0))
      else  
        Cs = Cs0
      end if
      vis_t(c) = density &
              * (lf*lf)           &          ! delta^2 
              * (Cs*Cs)           &          ! Cs^2   
              * shear(c) 
    end do

  else if(turbulence_model_variant == DYNAMIC) then
    if(buoyancy == YES) then  
      do c = 1, grid % n_cells
        lf = grid % vol(c)**ONE_THIRD  
        vis_t(c) = density            &
                * (lf*lf)             &          ! delta^2 
                * c_dyn(c)            &          ! c_dynamic   
                * shear(c) 
      end do
    else
      do c = 1, grid % n_cells 
        vis_t(c) = density                &
                * (lf*lf)                 &          ! delta^2 
                * c_dyn(c)                &          ! c_dynamic   
                * sqrt(shear(c)*shear(c)  &
                + 2.5*(grav_x*t_x(c) + grav_y*t_y(c) + grav_z*t_z(c)))  
      end do
    end if     

  else if(turbulence_model_variant == WALE) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD    
      vis_t(c) = density           &
              * (lf*lf)            &          ! delta^2 
              * (0.5*0.5)          &          ! Cs^2   
              * wale_v(c) 
    end do
  end if

  if(buoyancy == YES) then
    do c = 1, grid % n_cells
      Nc2 = -(  grav_x * t_x(c)   &
              + grav_y * t_y(c)   &
              + grav_z * t_z(c))  &
          / Tref
      Nc2 = max(0.0, Nc2) 
      vis_t(c) = vis_t(c) * sqrt(1.0 - min(2.5*Nc2/(shear(c)**2), 1.0))
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

      s_tot = sqrt(grid % sx(s)*grid % sx(s) +  &
                  grid % sy(s)*grid % sy(s) +  &
                  grid % sz(s)*grid % sz(s))

      nx = grid % sx(s) / s_tot 
      ny = grid % sy(s) / s_tot 
      nz = grid % sz(s) / s_tot 

      if(Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2)==WALLFL) then

        u_tot = sqrt(  u % n(c1) * u % n(c1)     &
                     + v % n(c1) * v % n(c1)     & 
                     + w % n(c1) * w % n(c1)  )

        u_nor = (   u % n(c1) * nx  &
                  + v % n(c1) * ny  &
                  + w % n(c1) * nz )   

        if( abs(u_tot) > abs(u_nor) ) then
          u_tan = sqrt(u_tot * u_tot - u_nor * u_nor)
        else
          u_tan = TINY 
        end if

        a_pow = 8.3
        b_pow = 1.0/7.0
        nu = viscosity/density
        dely = grid % wall_dist(c1)

        ! Calculate UtauL
        UtauL = ( u_tan/a_pow * (nu/dely)**b_pow )                     &
          ** (1.0/(1.0+b_pow))

        ! Calculate tau_wall 
        tau_wall(c1) = viscosity * u_tan / dely 
 
        ! Calculate y+
        y_plus  = dely*UtauL/nu
        if(y_plus  >=  11.81) then
          ! This one is effective viscosity
          VISwall(c1) = density*UtauL*UtauL*dely/abs(u_tan) 
        else 
          VISwall(c1) = viscosity + fw(s)*vis_t(c1)+(1.0-fw(s))*vis_t(c2)
        endif
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Exchange(grid, vis_t)
  call Exchange(grid, VISwall)

  end subroutine