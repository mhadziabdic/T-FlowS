!==============================================================================!
  subroutine Refine_Grid
!------------------------------------------------------------------------------!
!   Mark the region of the domain for local refinement and refine the grid!    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, lev, reg, n1, n2, n3, n4, n5, n6, n7, n8
  real                 :: x1, y1, z1, x8, y8, z8, x0, y0, z0
!==============================================================================!

  ! Set no cell for refinement, intially
  CelMar = 0

  do lev = 1,n_refine_levels

    do reg = 1,n_refined_regions(lev) 

      x1=refined_regions(lev,reg,1)
      y1=refined_regions(lev,reg,2)
      z1=refined_regions(lev,reg,3)
      x8=refined_regions(lev,reg,4)
      y8=refined_regions(lev,reg,5)
      z8=refined_regions(lev,reg,6)

      do c=1,NC
        n1 = grid % cells(c) % n(1)
        n2 = grid % cells(c) % n(2)
        n3 = grid % cells(c) % n(3)
        n4 = grid % cells(c) % n(4)
        n5 = grid % cells(c) % n(5)
        n6 = grid % cells(c) % n(6)
        n7 = grid % cells(c) % n(7)
        n8 = grid % cells(c) % n(8)

        x0=1.25e-1*(grid % nodes(n1) % x + grid % nodes(n2) % x +   &
                    grid % nodes(n3) % x + grid % nodes(n4) % x +   &
                    grid % nodes(n5) % x + grid % nodes(n6) % x +   &
                    grid % nodes(n7) % x + grid % nodes(n8) % x)
        y0=1.25e-1*(grid % nodes(n1) % y + grid % nodes(n2) % y +   &
                    grid % nodes(n3) % y + grid % nodes(n4) % y +   &
                    grid % nodes(n5) % y + grid % nodes(n6) % y +   &
                    grid % nodes(n7) % y + grid % nodes(n8) % y)
        z0=1.25e-1*(grid % nodes(n1) % z + grid % nodes(n2) % z +   &
                    grid % nodes(n3) % z + grid % nodes(n4) % z +   &
                    grid % nodes(n5) % z + grid % nodes(n6) % z +   &
                    grid % nodes(n7) % z + grid % nodes(n8) % z)

        if(refined_regions(lev,reg,0) == ELIPSOID) then
          if(  ( ((x1-x0)/x8)**2 +                                  &
                 ((y1-y0)/y8)**2 +                                  &
                 ((z1-z0)/z8)**2)  < 1.0 ) then
            CelMar(c) = -1
          end if
        else if(refined_regions(lev,reg,0) == RECTANGLE) then 
          if( (x1  < x0) .and. (x0  < x8) .and.                     &
              (y1  < y0) .and. (y0  < y8) .and.                     &
              (z1  < z0) .and. (z0  < z8) ) then
            CelMar(c) = -1
          endif
        else if(refined_regions(lev,reg,0) == PLANE) then 
          if( (x0-x1)*x8+(y0-y1)*y8+(z0-z1)*z8   >  0.0 ) then
            CelMar(c) = -1
          endif
        end if 
      end do   ! cells

    end do   ! reg

    call Refine_Marked_Cells(lev)

    CelMar = 0

  end do  ! lev

  end subroutine Refine_Grid
