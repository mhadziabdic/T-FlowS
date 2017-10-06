!======================================================================!
  subroutine Mark
!----------------------------------------------------------------------!
!   Mark the region of the domain for local refinement.                !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, lev, reg, n1, n2, n3, n4, n5, n6, n7, n8
  real    :: x1, y1, z1, x8, y8, z8, x0, y0, z0
!======================================================================!

  do c=-MAXB,MAXN
    CelMar(c) = 0
  end do 

  do lev = 1,NRL
    do reg = 1,n_refined_regions(lev) 
      x1=refined_regions(lev,reg,1)
      y1=refined_regions(lev,reg,2)
      z1=refined_regions(lev,reg,3)
      x8=refined_regions(lev,reg,4)
      y8=refined_regions(lev,reg,5)
      z8=refined_regions(lev,reg,6)

      do c=1,NC
        n1=CellN(c,1)
        n2=CellN(c,2)
        n3=CellN(c,3)
        n4=CellN(c,4)
        n5=CellN(c,5)
        n6=CellN(c,6)
        n7=CellN(c,7)
        n8=CellN(c,8)

        x0=1.25e-1*(x_node(n1)+x_node(n2)+x_node(n3)+x_node(n4)+  &
                    x_node(n5)+x_node(n6)+x_node(n7)+x_node(n8))
        y0=1.25e-1*(y_node(n1)+y_node(n2)+y_node(n3)+y_node(n4)+  &
                    y_node(n5)+y_node(n6)+y_node(n7)+y_node(n8))
        z0=1.25e-1*(z_node(n1)+z_node(n2)+z_node(n3)+z_node(n4)+  &
                    z_node(n5)+z_node(n6)+z_node(n7)+z_node(n8))

!->>>   write(*,*) x0, y0, z0

        if(refined_regions(lev,reg,0) == ELIPSO) then
          if(  ( ((x1-x0)/x8)**2 +                                  &
                 ((y1-y0)/y8)**2 +                                  &
                 ((z1-z0)/z8)**2)  < 1.0 ) then
            CelMar(c) = -1
          end if
        else if(refined_regions(lev,reg,0) == RECTAN) then 
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
      end do   ! =-> cells

    end do   ! =-> reg
    call refine(lev)

    do c=-MAXB,MAXN
       CelMar(c) = 0
    enddo 

  end do  ! =-> lev

  end subroutine Mark
