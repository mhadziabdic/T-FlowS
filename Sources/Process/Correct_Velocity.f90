!==============================================================================!
  real function Correct_Velocity(grid)
!------------------------------------------------------------------------------!
!   Corrects the velocities, and mass fluxes on the cell faces.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use Grid_Mod
  use Bulk_Mod
  use Info_Mod
  use Parameters_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer   :: c, c1, c2, s, m
  real      :: cfl_max(256), pe_max(256)
  real      :: cfl_t, pe_t
  real      :: Pdrop, FluxM
!==============================================================================!

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !- - - - - - - - - - - - - - - - - - - -  !
  !   For SOLIDs, px, py and pz are zero    !
  !   so this loop will not correct SOLID   !
  !   velocities.                           !
  !-----------------------------------------!
  if(ALGOR == FRACT) then
    do c = 1, grid % n_cells
      U % n(c) = U % n(c) - p % x(c) * grid % vol(c) / A % sav(c)
      V % n(c) = V % n(c) - p % y(c) * grid % vol(c) / A % sav(c)
      W % n(c) = W % n(c) - p % z(c) * grid % vol(c) / A % sav(c)
    end do 
  else ! algorythm is SIMPLE
    do c = 1, grid % n_cells
      U % n(c) = U % n(c) - p % x(c) * grid % vol(c) / A % sav(c)
      V % n(c) = V % n(c) - p % y(c) * grid % vol(c) / A % sav(c)
      W % n(c) = W % n(c) - p % z(c) * grid % vol(c) / A % sav(c)
    end do 
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  < 0) then
      if( (TypeBC(c2) == PRESSURE) ) then
        U % n(c2) = U % n(c1) 
        V % n(c2) = V % n(c1) 
        W % n(c2) = W % n(c1) 
      end if
    end if
  end do 

  !-------------------------------------------------------------------!
  !   Look at the following equation and you will understand why      !
  !   is the matrix for pressure corrections in SIMPLE algorythm      !
  !   formed from the coefficients of the velocity matrix.            !
  !   Moreover, it should also be clear that pressure correction      !
  !   matrix must be formed from underrelaxed velocity coefficients   !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
  !   Note that for FLUID-SOLID interaction, FLUX correctin is zero   !
  !   because A % val(A % pos(1,s)) is also zero.                     !  
  !   What will happen with parallel version ... only god knows.      !
  !-------------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
      if(c2  > 0) then
        Flux(s)=Flux(s)+(PP % n(c2) - PP % n(c1))*A % val(A % pos(1,s))
      else 
        Flux(s)=Flux(s)+(PP % n(c2) - PP % n(c1))*A % bou(c2)
      endif
    end if             !                                          !
  end do               !<---------- this is correction ---------->!

  !-------------------------------------!
  !    Calculate the max mass error     !
  !   with the new (corrected) fluxes   !
  !-------------------------------------!
  do c = 1, grid % n_cells
    b(c) = 0.0 
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
      b(c1)=b(c1)-Flux(s)
      if(c2  > 0) b(c2)=b(c2)+Flux(s)
    else
      b(c1) = b(c1)-Flux(s)
    end if
  end do

  do c = 1, grid % n_cells
    b(c) = b(c) / (grid % vol(c) * DENc(material(c)))
  end do

  errmax=0.0
  do c = 1, grid % n_cells
    errmax=max(errmax, abs(b(c)))
  end do
  call glomax(errmax)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  m = 1
  cfl_max(m) = 0.0
  pe_max(m)  = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if( (material(c1) .eq. m) .or. (material(c2) .eq. m) ) then
      if(c2  > 0 .or. c2  < 0.and.TypeBC(c2) == BUFFER) then
        cfl_t = abs( dt * Flux(s) /                &
                     ( Scoef(s) *                  &
                     (  grid % dx(s)*grid % dx(s)  &
                      + grid % dy(s)*grid % dy(s)  &
                      + grid % dz(s)*grid % dz(s)) ) )
        pe_t  = abs( Flux(s) / Scoef(s) / (VISc+TINY) )
        cfl_max(m) = max( cfl_max(m), cfl_t ) 
        pe_max(m)  = max( pe_max(m),  pe_t  ) 
      end if
    end if
  end do
  call glomax(cfl_max(m))
  call glomax(pe_max(m))

  call Info_Mod_Iter_Fill_At(1, 2, 'dum', -1, errmax)
  call Info_Mod_Bulk_Fill(cfl_max(m),          &
                          pe_max(m),           &                            
                          bulk(m) % flux_x,    &
                          bulk(m) % flux_y,    &
                          bulk(m) % flux_z,    &
                          bulk(m) % p_drop_x,  &
                          bulk(m) % p_drop_y,  &
                          bulk(m) % p_drop_z)

  Correct_Velocity = errmax ! /(velmax+TINY)

  end function
