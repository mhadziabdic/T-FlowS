!==============================================================================!
  subroutine Source_Kin_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
  use Work_Mod, only: kin_x => r_cell_01,  &
                      kin_y => r_cell_02,  &
                      kin_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, c1, c2, s 
  real              :: UtotSq, Unor, UnorSq, Utan
!==============================================================================!

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  if(turbulence_model_variant == HIGH_RE) then
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)
 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) /= BUFFER) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL) then

          ! Compute tangential velocity component
          UtotSq = u % n(c1) * u % n(c1) &
                 + v % n(c1) * v % n(c1) &
                 + w % n(c1) * w % n(c1)
          Unor = ( u % n(c1) * grid % sx(s)     &
                 + v % n(c1) * grid % sy(s)     &
                 + w % n(c1) * grid % sz(s) )   &
                 / sqrt(  grid % sx(s)*grid % sx(s)    &
                        + grid % sy(s)*grid % sy(s)    &
                        + grid % sz(s)*grid % sz(s))
          UnorSq = Unor*Unor

          if(UtotSq  >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if

          ! Compute nondimensional wall distance and wall-shear stress
          if(ROUGH == NO) then
            Ynd(c1) = sqrt(kin % n(c1)) * Cmu25 * grid % wall_dist(c1)  &
                    / viscosity
            tau_wall(c1) = abs(density                   &
                        * kappa * sqrt(kin % n(c1)) * Cmu25 * Utan &
                        / (log(Elog*Ynd(c1))))  

            ! Compute production in the first wall cell 
            p_kin(c1) = tau_wall(c1) * Cmu25 * sqrt(kin % n(c1)) &
                   / (kappa*grid % wall_dist(c1))

          else if(ROUGH==YES) then
            Ynd(c1) = sqrt(kin % n(c1)) * Cmu25 * (grid % wall_dist(c1)+Zo)  &
                    / viscosity
            tau_wall(c1) = abs(density                     &
                        * kappa * sqrt(kin % n(c1)) * Cmu25 * Utan   &
                        / (log((grid % wall_dist(c1)+Zo)/Zo)))  

            p_kin(c1) = tau_wall(c1) * Cmu25 * sqrt(kin % n(c1)) &
                   / (kappa*(grid % wall_dist(c1)+Zo))
            kin % n(c2) = tau_wall(c1)/0.09**0.5
          end if  

          ! Filling up the source term
          b(c1) = b(c1) + p_kin(c1) * grid % vol(c1) 
        end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do

    !-----------------------------------------!
    !   Compute the sources in the interior   !
    !-----------------------------------------!
    do c=1,grid % n_cells

      ! grid % cell_near_wall ensures not to put p_kin twice into the near wall cells
      if(.not. grid % cell_near_wall(c)) then

        ! Production:
        p_kin(c)= vis_t(c) * shear(c)*shear(c)
        b(c) = b(c) + p_kin(c) * grid % vol(c)
      end if

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c))                          &
                          + density * eps % n(c) / kin % n(c)  &
                          * grid % vol(c)
    end do
  end if    ! end if mode = wf 

  !--------------------------------------------------------!
  !   Jones-Launder model and Launder-Sharma + Yap model   !
  !--------------------------------------------------------!
  if(turbulence_model_variant == LOW_RE) then
    do c = 1, grid % n_cells

      ! Production:
      p_kin(c)= vis_t(c) * shear(c) * shear(c)
      b(c) = b(c) + p_kin(c) * grid % vol(c)

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c)) + density*eps%n(c)/(kin%n(c)+TINY)*grid % vol(c)

      ! Preparation of kin for the boundary condition. kin variable is temporaraly borrowed.
      kin % n(c) = sqrt(kin % n(c))
    end do

    call Grad_Mod_For_Phi(grid, kin % n, 1, kin_x, .true.)  ! dk/dx
    call Grad_Mod_For_Phi(grid, kin % n, 2, kin_y, .true.)  ! dk/dy
    call Grad_Mod_For_Phi(grid, kin % n, 3, kin_z, .true.)  ! dk/dz

    do c = 1, grid % n_cells

      ! Turning kin back to its real value
      kin % n(c) = kin % n(c) * kin % n(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))                &
                          + 2.0 * viscosity*(  kin_x(c)**2     &
                                             + kin_y(c)**2     &
                                             + kin_z(c)**2 )   &
                           * grid % vol(c) / (kin % n(c) + TINY)          
    end do
  end if

  if(turbulence_model == HYBRID_PITM) then
    do c = 1, grid % n_cells

      ! Production:
      p_kin(c)= vis_t(c) * shear(c) * shear(c)
      b(c) = b(c) + p_kin(c) * grid % vol(c)

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c))      &
                          + density                  &
                          * eps % n(c) / kin % n(c)  &
                          * grid % vol(c)
    end do
  end if

  call Comm_Mod_Exchange(grid, kin % n)

  end subroutine
