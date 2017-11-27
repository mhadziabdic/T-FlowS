!==============================================================================!
  subroutine Compute_Fluxes(grid)
!------------------------------------------------------------------------------!
!   Compute mass fluxes through whole domain.                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, m
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
!==============================================================================!

  do m = 1, grid % n_materials

    bulk(m) % flux_x = 0.0
    bulk(m) % flux_y = 0.0
    bulk(m) % flux_z = 0.0

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 > 0) then
        if( (material(c1) == m) .and.  &
            (material(c1) == material(c2)) ) then
          xc1=grid % xc(c1) 
          yc1=grid % yc(c1) 
          zc1=grid % zc(c1) 
          xc2=grid % xc(c1) + grid % dx(s) 
          yc2=grid % yc(c1) + grid % dy(s) 
          zc2=grid % zc(c1) + grid % dz(s)

          if((xc1 <= xp(m)).and.(xc2 > xp(m)))  &
            bulk(m) % flux_x = bulk(m) % flux_x + Flux(s)

          if((yc1 <= yp(m)).and.(yc2 > yp(m)))  &
            bulk(m) % flux_y = bulk(m) % flux_y + Flux(s)

          if((zc1 <= zp(m)).and.(zc2 > zp(m)))  &
            bulk(m) % flux_z = bulk(m) % flux_z + Flux(s)

          if((xc2 < xp(m)).and.(xc1 >= xp(m)))  &
            bulk(m) % flux_x = bulk(m) % flux_x - Flux(s)

          if((yc2 < yp(m)).and.(yc1 >= yp(m)))  &
            bulk(m) % flux_y = bulk(m) % flux_y - Flux(s)

          if((zc2 < zp(m)).and.(zc1 >= zp(m)))  &
            bulk(m) % flux_z = bulk(m) % flux_z - Flux(s)

        end if ! material 1&2
      else if(c2 < 0.and.TypeBC(c2) == BUFFER) then
        if( (material(c1)==m) .and. (material(c1) == material(c2)) ) then
          xc1=grid % xc(c1) 
          yc1=grid % yc(c1) 
          zc1=grid % zc(c1) 
          xc2=grid % xc(c1) + grid % dx(s) 
          yc2=grid % yc(c1) + grid % dy(s) 
          zc2=grid % zc(c1) + grid % dz(s)

          if((xc1 <= xp(m)).and.(xc2 > xp(m)))  &
            bulk(m) % flux_x = bulk(m) % flux_x + .5*Flux(s)

          if((yc1 <= yp(m)).and.(yc2 > yp(m)))  &
            bulk(m) % flux_y = bulk(m) % flux_y + .5*Flux(s)

          if((zc1 <= zp(m)).and.(zc2 > zp(m)))  &
            bulk(m) % flux_z = bulk(m) % flux_z + .5*Flux(s)

          if((xc2 < xp(m)).and.(xc1 >= xp(m)))  &
            bulk(m) % flux_x = bulk(m) % flux_x - .5*Flux(s)

          if((yc2 < yp(m)).and.(yc1 >= yp(m)))  &
            bulk(m) % flux_y = bulk(m) % flux_y - .5*Flux(s)

          if((zc2 < zp(m)).and.(zc1 >= zp(m)))  &
            bulk(m) % flux_z = bulk(m) % flux_z - .5*Flux(s)

        end if ! material 1&2
      end if   ! c2 > 0
    end do

    call glosum(bulk(m) % flux_x)
    call glosum(bulk(m) % flux_y)
    call glosum(bulk(m) % flux_z)

    bulk(m) % u = bulk(m) % flux_x / (bulk(m) % area_x + TINY)
    bulk(m) % v = bulk(m) % flux_y / (bulk(m) % area_y + TINY)
    bulk(m) % w = bulk(m) % flux_z / (bulk(m) % area_z + TINY)

  end do ! m

  end subroutine
