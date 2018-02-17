!==============================================================================!
  subroutine Source_Kin_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Parameters_Mod
  use Work_Mod, only: kin_x => r_cell_01,  &
                      kin_y => r_cell_02,  &
                      kin_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s 
  real    :: UtotSq, Unor, UnorSq, Utan
!==============================================================================!

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  if(MODE == HIGH_RE) then
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)
 
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

          ! Compute tangential velocity component
          UtotSq = U % n(c1) * U % n(c1) &
                 + V % n(c1) * V % n(c1) &
                 + W % n(c1) * W % n(c1)
          Unor = ( U % n(c1) * grid % sx(s)     &
                 + V % n(c1) * grid % sy(s)     &
                 + W % n(c1) * grid % sz(s) )   &
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
            Ynd(c1) = sqrt(kin % n(c1)) * Cmu25 * WallDs(c1)  &
                    / VISc
            TauWall(c1) = abs(DENc(material(c1))                   &
                        * kappa * sqrt(kin % n(c1)) * Cmu25 * Utan &
                        / (log(Elog*Ynd(c1))))  

            ! Compute production in the first wall cell 
            Pk(c1) = TauWall(c1) * Cmu25 * sqrt(kin % n(c1)) &
                   / (kappa*WallDs(c1))

          else if(ROUGH==YES) then
            Ynd(c1) = sqrt(kin % n(c1)) * Cmu25 * (WallDs(c1)+Zo)  &
                    / VISc
            TauWall(c1) = abs(DENc(material(c1))                     &
                        * kappa * sqrt(kin % n(c1)) * Cmu25 * Utan   &
                        / (log((WallDs(c1)+Zo)/Zo)))  

            Pk(c1) = TauWall(c1) * Cmu25 * sqrt(kin % n(c1)) &
                   / (kappa*(WallDs(c1)+Zo))
            kin % n(c2) = TauWall(c1)/0.09**0.5
          end if  

          ! Filling up the source term
          b(c1) = b(c1) + Pk(c1) * grid % vol(c1) 
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do

    !-----------------------------------------!
    !   Compute the sources in the interior   !
    !-----------------------------------------!
    do c=1,grid % n_cells

      ! IsNearWall ensures not to put Pk twice into the near wall cells
      if(.NOT. IsNearWall(c)) then

        ! Production:
        Pk(c)= VISt(c) * Shear(c)*Shear(c)
        b(c) = b(c) + Pk(c) * grid % vol(c)
      end if

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c))                          &
                          + DENc(material(c)) * eps % n(c) / kin % n(c)  &
                          * grid % vol(c)
    end do
  end if    ! end if mode = wf 

  !--------------------------------------------------------!
  !   Jones-Launder model and Launder-Sharma + Yap model   !
  !--------------------------------------------------------!
  if(MODE == LOW_RE) then
    do c = 1, grid % n_cells

      ! Production:
      Pk(c)= VISt(c) * Shear(c) * Shear(c)
      b(c) = b(c) + Pk(c) * grid % vol(c)

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c)) + &
      DENc(material(c))*eps%n(c)/(kin%n(c)+TINY)*grid % vol(c)

      if (BUOY == YES) then
        b(c) = b(c) + Pbuoy(c) * grid % vol(c)
      end if

      ! Preparation of kin for the boundary condition. kin variable is temporaraly borrowed.
      kin % n(c) = sqrt(kin % n(c))
    end do

    call GraPhi(kin % n, 1, kin_x, .TRUE.)  ! dK/dx
    call GraPhi(kin % n, 2, kin_y, .TRUE.)  ! dK/dy
    call GraPhi(kin % n, 3, kin_z, .TRUE.)  ! dK/dz

    do c = 1, grid % n_cells

      ! Turning kin back to its real value
      kin % n(c) = kin % n(c) * kin % n(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))           &
                          + 2.0 * VISc*(  kin_x(c)**2     &
                                        + kin_y(c)**2     &
                                        + kin_z(c)**2 )   &
                           * grid % vol(c) / (kin % n(c) + TINY)          
    end do
  end if

  if(SIMULA == HYB_PITM) then
    do c = 1, grid % n_cells

      ! Production:
      Pk(c)= VISt(c) * Shear(c) * Shear(c)
      b(c) = b(c) + Pk(c) * grid % vol(c)

      ! Dissipation:
      A % val(A % dia(c)) = A % val(A % dia(c))      &
                          + DENc(material(c))        &
                          * eps % n(c) / kin % n(c)  &
                          * grid % vol(c)
    end do
  end if

  call Exchange(grid, kin % n)

  end subroutine
