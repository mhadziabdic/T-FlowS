!==============================================================================!
  subroutine Compute_Grid_Geometry(rrun)
!------------------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(in) :: rrun     
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume        
  real :: Distance       
  real :: Distance_Squared       
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, m, n, s, c_1, c_2
  integer :: wall_mark
  real    :: local_x_node(4), local_y_node(4), local_z_node(4)
  real    :: x_cell_tmp, y_cell_tmp, z_cell_tmp    
  real    :: xs2, ys2, zs2
  real    :: dsc1, dsc2          !  for the interpolation factors
  real    :: t, tot_surf , maxdis
  real    :: xc1, yc1, zc1, xc2, yc2, zc2 
  real    :: x_min, x_max, y_min, y_max, z_min, z_max
  integer :: f4n(6,4)
  integer :: f3n(4,3)
!==============================================================================!
!
!                                n3 
!                 +---------------!---------------+
!                /|              /|              /|
!               / |             / |             / |
!              /  |          n2/  |            /  |
!             +---------------!---------------+   |
!             |   |           |   |           |   |
!             |   |     o---- | s-------o     |   |      
!             |   +---- c1 ---|   !---- c2 ---|   +       
!             |  /            |  /n4          |  /
!             | /             | /             | /
!             |/              |/              |/
!             +---------------!---------------+
!                            n1
!
!   Notes:
! 
!     ! side s is oriented from cell center c1 to cell center c2     
!     ! c2 is greater then c1 inside the domain or smaller then 0
!       on the boundary
!     ! nodes are denoted with n1 - n4
!
!            c3           
!             \  4-----3
!              \/ \  . |
!              /   \  +---c2
!             /  .  \  |
!            / .     \ |
!           /.        \|
!          1-----------2
!                   |
!                   c1
!
!                                n3 
!                 +---------------!-------+         
!                /|            n2/|      /|
!               / |             !-------+ |
!              /  |            /|s|  c2 | |
!             +---------------+ | !-----| +
!             |   |           | |/n4    |/
!             |   |     c1    | !-------+
!             |   +-----------|n1 +
!             |  /            |  /
!             | /             | /
!             |/              |/ 
!             +---------------+
!                            n1
! 
!------------------------------------------------------------------------------!
!   Generaly:
!
!   the equation of plane reads: A*x_node + B*y_node + C*z_node + D = 0
!
!   and the equation of line:  x_node = x0 + t*rx
!                              y_node = y0 + t*ry
!                              z_node = z0 + t*rz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   In our case:
!
!     line is a connection between the two cell centers:
!
!     x_node = xc(c1) + t*(xc(c2)-xc(c1)) = xc(c1) + t*rx
!     y_node = yc(c1) + t*(yc(c2)-yc(c1)) = yc(c1) + t*ry
!     z_node = zc(c1) + t*(zc(c2)-zc(c1)) = zc(c1) + t*rz
!    
!
!     plane is a cell face: 
!
!     Sx * x_node + Sy * y_node + Sz * z_node = Sx * xsp(s) + Sy * ysp(s) + Sz * zsp(s)
!  
!     and the intersection is at:
!  
!         Sx*(xsp(s)-xc(c1)) + Sy*(ysp(s)-yc(c1) + Sz*(zsp(s)-zc(c1)) 
!     t = -----------------------------------------------------------
!                           rx*Sx + ry*Sy + rz*Sz
!  
!------------------------------------------------------------------------------!
  data    f4n / 1, 1, 2, 4, 3, 5,                                   &
                2, 5, 6, 8, 7, 7,                                   &
                4, 6, 8, 7, 5, 8,                                   &
                3, 2, 4, 3, 1, 6  /

  data    f3n / 1,  1,  2,  3,                                      &
                2,  4,  4,  4,                                      &
                3,  2,  3,  1 /

    ! Without the following six lines, this procedure works for any grid
    do c=1,NC
      CellN(c,0)=8
    end do
    do s=1,NS
      SideN(s,0)=4 
    end do

    !-----------------------------------------!
    !   Calculate the cell centers            !
    !-----------------------------------------!
    !   => depends on: x_node,y_node,z_node   ! 
    !   <= gives:      xc,yc,zc c>0           !
    !-----------------------------------------!
    do c=1,NC
      xc(c)=0.0
      yc(c)=0.0
      zc(c)=0.0
      do n=1,CellN(c,0)
        xc(c) = xc(c) + x_node(CellN(c,n))/(1.0*CellN(c,0))
        yc(c) = yc(c) + y_node(CellN(c,n))/(1.0*CellN(c,0))
        zc(c) = zc(c) + z_node(CellN(c,n))/(1.0*CellN(c,0))
      end do
    end do

    !-----------------------------------------!
    !   Calculate delta                       !
    !-----------------------------------------!
    !   => depends on: x_node,y_node,z_node   ! 
    !   <= gives:      delta                  !
    !-----------------------------------------!
    do c=1,NC
      delta(c)=0.0
      x_min = +HUGE   
      y_min = +HUGE  
      z_min = +HUGE  
      x_max = -HUGE  
      y_max = -HUGE  
      z_max = -HUGE  
      do n=1,CellN(c,0)
        x_min = min(x_min, x_node(CellN(c,n)))
        y_min = min(y_min, y_node(CellN(c,n)))
        z_min = min(z_min, z_node(CellN(c,n)))
        x_max = max(x_max, x_node(CellN(c,n)))
        y_max = max(y_max, y_node(CellN(c,n)))
        z_max = max(z_max, z_node(CellN(c,n)))
      end do
      delta(c) = x_max-x_min
      delta(c) = max(delta(c), (y_max-y_min))
      delta(c) = max(delta(c), (z_max-z_min))
    end do

    !-----------------------------------------------------!
    !   Calculate:                                        ! 
    !      components of cell sides, cell side centers.   !
    !-----------------------------------------------------!
    !   => depends on: x_node,y_node,z_node               ! 
    !   <= gives:      Sx,Sy,Sz,xsp,yzp,zsp               !
    !-----------------------------------------------------!
    do s=1,NS
      do n=1,SideN(s,0)    ! for quadrilateral an triangular faces
        local_x_node(n)=x_node(SideN(s,n))
        local_y_node(n)=y_node(SideN(s,n))
        local_z_node(n)=z_node(SideN(s,n))
      end do                       

      ! Cell side components
      if( SideN(s,0)  ==  4 ) then
        Sx(s)= 0.5 * (   (local_y_node(2)-local_y_node(1))  &
                       * (local_z_node(2)+local_z_node(1))  &
                       + (local_y_node(3)-local_y_node(2))  &
                       * (local_z_node(2)+local_z_node(3))  &
                       + (local_y_node(4)-local_y_node(3))  &
                       * (local_z_node(3)+local_z_node(4))  &
                       + (local_y_node(1)-local_y_node(4))  &
                       * (local_z_node(4)+local_z_node(1)) )
        Sy(s)= 0.5 * (   (local_z_node(2)-local_z_node(1))  &
                       * (local_x_node(2)+local_x_node(1))  &
                       + (local_z_node(3)-local_z_node(2))  &
                       * (local_x_node(2)+local_x_node(3))  &
                       + (local_z_node(4)-local_z_node(3))  &
                       * (local_x_node(3)+local_x_node(4))  &
                       + (local_z_node(1)-local_z_node(4))  &
                       * (local_x_node(4)+local_x_node(1)) )
        Sz(s)= 0.5 * (   (local_x_node(2)-local_x_node(1))  & 
                       * (local_y_node(2)+local_y_node(1))  & 
                       + (local_x_node(3)-local_x_node(2))  & 
                       * (local_y_node(2)+local_y_node(3))  &
                       + (local_x_node(4)-local_x_node(3))  & 
                       * (local_y_node(3)+local_y_node(4))  &
                       + (local_x_node(1)-local_x_node(4))  & 
                       * (local_y_node(4)+local_y_node(1)) )
      else if( SideN(s,0)  ==  3 ) then 
        Sx(s)= 0.5 * (   (local_y_node(2)-local_y_node(1))  &
                       * (local_z_node(2)+local_z_node(1))  & 
                       + (local_y_node(3)-local_y_node(2))  &
                       * (local_z_node(2)+local_z_node(3))  &
                       + (local_y_node(1)-local_y_node(3))  &
                       * (local_z_node(3)+local_z_node(1)) )
        Sy(s)= 0.5 * (   (local_z_node(2)-local_z_node(1))  &
                       * (local_x_node(2)+local_x_node(1))  &
                       + (local_z_node(3)-local_z_node(2))  &
                       * (local_x_node(2)+local_x_node(3))  & 
                       + (local_z_node(1)-local_z_node(3))  &
                       * (local_x_node(3)+local_x_node(1)) )
        Sz(s)= 0.5 * (   (local_x_node(2)-local_x_node(1))  &
                       * (local_y_node(2)+local_y_node(1))  &
                       + (local_x_node(3)-local_x_node(2))  & 
                       * (local_y_node(2)+local_y_node(3))  & 
                       + (local_x_node(1)-local_x_node(3))  & 
                       * (local_y_node(3)+local_y_node(1)) )
      else
        write(*,*) 'Compute_Grid_Geometry: something horrible has happened !'
        stop
      end if

      ! Barycenters
      if(SideN(s,0) == 4) then  
        xsp(s) = (   local_x_node(1)+local_x_node(2)        &
                   + local_x_node(3)+local_x_node(4) )/4.0
        ysp(s) = (   local_y_node(1)+local_y_node(2)        &
                   + local_y_node(3)+local_y_node(4) )/4.0
        zsp(s) = (   local_z_node(1)+local_z_node(2)        &
                   + local_z_node(3)+local_z_node(4) )/4.0
      else if(SideN(s,0) == 3) then  
        xsp(s) = (local_x_node(1)+local_x_node(2)+local_x_node(3))/3.0
        ysp(s) = (local_y_node(1)+local_y_node(2)+local_y_node(3))/3.0
        zsp(s) = (local_z_node(1)+local_z_node(2)+local_z_node(3))/3.0
      end if 

    end do ! through sides

    !--------------------------------------!
    !   Calculate boundary cell centers    !
    !--------------------------------------!
    !   => depends on: xc,yc,zc,Sx,Sy,Sz   !
    !   <= gives:      xc,yc,zc for c<0    !   
    !--------------------------------------!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      tot_surf = sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))

      if(c2  < 0) then
        t = (   Sx(s)*(xsp(s)-xc(c1))                               &
              + Sy(s)*(ysp(s)-yc(c1))                               &
              + Sz(s)*(zsp(s)-zc(c1)) ) / tot_surf
        xc(c2) = xc(c1) + Sx(s)*t / tot_surf
        yc(c2) = yc(c1) + Sy(s)*t / tot_surf
        zc(c2) = zc(c1) + Sz(s)*t / tot_surf
      endif 
    end do ! through sides

    !---------------------------------------------!
    !   Find the sides on the periodic boundary   !
    !---------------------------------------------!
    !   => depends on: xc,yc,zc,Sx,Sy,Sz          !
    !   <= gives:      Dx,Dy,Dz                   !
    !---------------------------------------------!
    if(rrun) then
    NSsh = 0
    do s=1,NS

      ! Initialize
      Dx(s)=0.0
      Dy(s)=0.0
      Dz(s)=0.0

      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2   >  0) then

        ! Scalar product of the side with line c1-c2 is a good criterion
        if( (Sx(s) * (xc(c2)-xc(c1) )+                              &
             Sy(s) * (yc(c2)-yc(c1) )+                              &
             Sz(s) * (zc(c2)-zc(c1) ))  < 0.0 ) then

          NSsh = NSsh + 2
 
          ! Find the coordinates of ...
          m=SideCc(s,2)

          if(SideN(s,0) == 4) then   

            ! Coordinates of the shadow face
            xs2=.25*(x_node(CellN(c2,f4n(m,1)))+x_node(CellN(c2,f4n(m,2)))+   &
                     x_node(CellN(c2,f4n(m,3)))+x_node(CellN(c2,f4n(m,4))) )
            ys2=.25*(y_node(CellN(c2,f4n(m,1)))+y_node(CellN(c2,f4n(m,2)))+   &
                     y_node(CellN(c2,f4n(m,3)))+y_node(CellN(c2,f4n(m,4))) )
            zs2=.25*(z_node(CellN(c2,f4n(m,1)))+z_node(CellN(c2,f4n(m,2)))+   &
                     z_node(CellN(c2,f4n(m,3)))+z_node(CellN(c2,f4n(m,4))) )
 
            ! Add shadow faces
            SideN(NS+NSsh-1,0) = 4
            SideC(1,NS+NSsh-1) = c1 
            SideC(2,NS+NSsh-1) = -NbC-1
            SideN(NS+NSsh-1,1) = SideN(s,1)
            SideN(NS+NSsh-1,2) = SideN(s,2)
            SideN(NS+NSsh-1,3) = SideN(s,3)
            SideN(NS+NSsh-1,4) = SideN(s,4)
            Sx(NS+NSsh-1) = Sx(s)
            Sy(NS+NSsh-1) = Sy(s)
            Sz(NS+NSsh-1) = Sz(s)
            xsp(NS+NSsh-1) = xsp(s)
            ysp(NS+NSsh-1) = ysp(s)
            zsp(NS+NSsh-1) = zsp(s)
            SideN(NS+NSsh,0) = 4
            SideC(1,NS+NSsh) = c2 
            SideC(2,NS+NSsh) = -NbC-1
            SideN(NS+NSsh,1) = CellN(c2,f4n(m,1)) 
            SideN(NS+NSsh,2) = CellN(c2,f4n(m,2))
            SideN(NS+NSsh,3) = CellN(c2,f4n(m,3))
            SideN(NS+NSsh,4) = CellN(c2,f4n(m,4))
            Sx(NS+NSsh) = Sx(s)
            Sy(NS+NSsh) = Sy(s)
            Sz(NS+NSsh) = Sz(s)
            xsp(NS+NSsh) = xs2
            ysp(NS+NSsh) = ys2
            zsp(NS+NSsh) = zs2
          else if(SideN(s,0) == 3) then  

            ! Coordinates of the shadow face
            xs2=.33333333 * (x_node(CellN(c2,f3n(m,1)))+                 &
                     x_node(CellN(c2,f3n(m,2)))+x_node(CellN(c2,f3n(m,3))) )
            ys2=.33333333 * (y_node(CellN(c2,f4n(m,1)))+                 &
                     y_node(CellN(c2,f3n(m,2)))+y_node(CellN(c2,f3n(m,3))) )
            zs2=.33333333 * (z_node(CellN(c2,f4n(m,1)))+                 &
                     z_node(CellN(c2,f3n(m,2)))+z_node(CellN(c2,f3n(m,3))) )

            ! Add shadow faces
            SideN(NS+NSsh-1,0) = 3
            SideC(1,NS+NSsh-1) = c1 
            SideC(2,NS+NSsh-1) = -NbC-1
            SideN(NS+NSsh-1,1) = SideN(s,1)
            SideN(NS+NSsh-1,2) = SideN(s,2)
            SideN(NS+NSsh-1,3) = SideN(s,3)
            Sx(NS+NSsh-1) = Sx(s)
            Sy(NS+NSsh-1) = Sy(s)
            Sz(NS+NSsh-1) = Sz(s)
            xsp(NS+NSsh-1) = xsp(s)
            ysp(NS+NSsh-1) = ysp(s)
            zsp(NS+NSsh-1) = zsp(s)
            SideN(NS+NSsh,0) = 3
            SideC(1,NS+NSsh) = c2 
            SideC(2,NS+NSsh) = -NbC-1
            SideN(NS+NSsh,1) = CellN(c2,f3n(m,1)) 
            SideN(NS+NSsh,2) = CellN(c2,f3n(m,2))
            SideN(NS+NSsh,3) = CellN(c2,f3n(m,3))
            Sx(NS+NSsh) = Sx(s)
            Sy(NS+NSsh) = Sy(s)
            Sz(NS+NSsh) = Sz(s)
            xsp(NS+NSsh) = xs2
            ysp(NS+NSsh) = ys2
            zsp(NS+NSsh) = zs2
          end if 

          Dx(s)=xsp(s)-xs2  !------------------------!
          Dy(s)=ysp(s)-ys2  ! later: xc2 = xc2 + Dx  !
          Dz(s)=zsp(s)-zs2  !------------------------!

        endif !  S*(c2-c1) < 0.0
      end if  !  c2 > 0
    end do    !  sides  
    write(*,*) 'Number of shadow faces: ', NSsh
    end if

  !----------------------------------!
  !   Calculate the cell volumes     !
  !----------------------------------!
  !   => depends on: xc,yc,zc,       !
  !                  Dx,Dy,Dz,       !
  !                  xsp, ysp, zsp   !
  !   <= gives:      volume          !
  !----------------------------------!
  if(rrun) then
  do c=1,NC
    volume(c)=0.0
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)   

    do n=1,SideN(s,0)      ! for quadrilateral an triangular faces
      local_x_node(n)=x_node(SideN(s,n))
      local_y_node(n)=y_node(SideN(s,n))
      local_z_node(n)=z_node(SideN(s,n))
    end do   

    ! First cell
    x_cell_tmp=xc(c1)
    y_cell_tmp=yc(c1)
    z_cell_tmp=zc(c1)
    dsc1=Distance(x_cell_tmp,y_cell_tmp,z_cell_tmp,xsp(s), ysp(s), zsp(s)) 
    volume(c1)=volume(c1) + Tet_Volume(xsp(s),ysp(s),zsp(s),            &
               local_x_node(1),local_y_node(1),local_z_node(1),     &
               local_x_node(2),local_y_node(2),local_z_node(2),     &
               x_cell_tmp,y_cell_tmp,z_cell_tmp)
    volume(c1)=volume(c1) + Tet_Volume(xsp(s),ysp(s),zsp(s),            &
               local_x_node(2),local_y_node(2),local_z_node(2),     &
               local_x_node(3),local_y_node(3),local_z_node(3),     &
               x_cell_tmp,y_cell_tmp,z_cell_tmp)
    if(SideN(s,0) == 4) then
      volume(c1)=volume(c1) + Tet_Volume(xsp(s),ysp(s),zsp(s),          &
                 local_x_node(3),local_y_node(3),local_z_node(3),   &
                 local_x_node(4),local_y_node(4),local_z_node(4),   &
                 x_cell_tmp,y_cell_tmp,z_cell_tmp)
      volume(c1)=volume(c1) + Tet_Volume(xsp(s),ysp(s),zsp(s),          &
                 local_x_node(4),local_y_node(4),local_z_node(4),   &
                 local_x_node(1),local_y_node(1),local_z_node(1),   &
                 x_cell_tmp,y_cell_tmp,z_cell_tmp)
    else if(SideN(s,0) == 3) then
      volume(c1)=volume(c1) + Tet_Volume(xsp(s),ysp(s),zsp(s),          &
                 local_x_node(3),local_y_node(3),local_z_node(3),   &
                 local_x_node(1),local_y_node(1),local_z_node(1),   &
                 x_cell_tmp,y_cell_tmp,z_cell_tmp)
    end if

    ! Second cell
    if(c2  > 0) then
      x_cell_tmp=xc(c2)+Dx(s)
      y_cell_tmp=yc(c2)+Dy(s)
      z_cell_tmp=zc(c2)+Dz(s)
      dsc2=Distance(x_cell_tmp,y_cell_tmp,z_cell_tmp,xsp(s), ysp(s), zsp(s)) 
      volume(c2)=volume(c2) -Tet_Volume(xsp(s),ysp(s),zsp(s),           &
                 local_x_node(1),local_y_node(1),local_z_node(1),   &
                 local_x_node(2),local_y_node(2),local_z_node(2),   &
                 x_cell_tmp,y_cell_tmp,z_cell_tmp)
      volume(c2)=volume(c2) -Tet_Volume(xsp(s),ysp(s),zsp(s),           &
                 local_x_node(2),local_y_node(2),local_z_node(2),   &
                 local_x_node(3),local_y_node(3),local_z_node(3),   &
                 x_cell_tmp,y_cell_tmp,z_cell_tmp)
      if(SideN(s,0) == 4) then
        volume(c2)=volume(c2) -Tet_Volume(xsp(s),ysp(s),zsp(s),         &
                   local_x_node(3),local_y_node(3),local_z_node(3), &
                   local_x_node(4),local_y_node(4),local_z_node(4), &
                   x_cell_tmp,y_cell_tmp,z_cell_tmp)
        volume(c2)=volume(c2) -Tet_Volume(xsp(s),ysp(s),zsp(s),         &
                   local_x_node(4),local_y_node(4),local_z_node(4), &
                   local_x_node(1),local_y_node(1),local_z_node(1), &
                   x_cell_tmp,y_cell_tmp,z_cell_tmp)
      else if(SideN(s,0) == 3) then
        volume(c2)=volume(c2) -Tet_Volume(xsp(s),ysp(s),zsp(s),         &
                   local_x_node(3),local_y_node(3),local_z_node(3), &
                   local_x_node(1),local_y_node(1),local_z_node(1), &
                   x_cell_tmp,y_cell_tmp,z_cell_tmp)
      end if  
    else        
      dsc2=0.0
    end if

  end do
  end if

  !----------------------------------------------------------!
  !   Calculate:                                             ! 
  !      distance from the cell center to the nearest wall   !
  !----------------------------------------------------------!
  !   => depends on: xc,yc,zc inside and on the boundary     !
  !   <= gives:      WallDs i                                !
  !----------------------------------------------------------!
  WallDs = HUGE 

  write(*,'(A43)')   '#------------------------------------------'
  if(rrun) then
    write(*,'(A43)') '# Computing the distance to the walls (2/2)'           
  else            
    write(*,'(A43)') '# Computing the distance to the walls (1/2)'           
  end if 
  write(*,'(A43)')   '#------------------------------------------'
  write(*,*) 'Insert the highest BC marker which represents the wall'
  write(*,*) 'for computing the distance to the wall (0 to skip)'
  write(*,*) 'If you skip this, smoothing will not work properly'
  read(*,*) wall_mark   

  if(wall_mark == 0) then
    WallDs = 1.0
    write(*,*) 'Distance to the wall set to 1 everywhere !'            
  else 
    do c1=1,NC 
      do s = WallFacFst, WallFacLst      ! 1,NS
        c_1 = SideC(1,s)
        c_2 = SideC(2,s)
        if(c_2 < 0) then
          if(BCmark(c_2) <= wall_mark) then
            WallDs(c1)=min(WallDs(c1), &
            Distance_Squared(xc(c1),yc(c1),zc(c1),xsp(s),ysp(s),zsp(s)))
          end if
!        else if(material(c_1) /= material(c_2) ) then
!          WallDs(c1)=min(WallDs(c1), &
!          Distance_Squared(xc(c1),yc(c1),zc(c1),xsp(s),ysp(s),zsp(s)))
        end if 
      end do
    end do

    do c=1,NC
      WallDs(c)=sqrt(WallDs(c))
    end do

    write(*,*) 'Maximal distance to the wall: ', maxval(WallDs(1:NC))
    write(*,*) 'Minimal distance to the wall: ', minval(WallDs(1:NC))
  end if

  do n=1,NN
    walln(n)=HUGE
  end do

  do c=1,NC
    do n=1,CellN(c,0)
      walln(CellN(c,n))=min(WallDs(c),walln(CellN(c,n)))
    end do
  end do

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0 .or. material(c1) /= material(c2)) then 
      do n=1,SideN(s,0)    ! for quadrilateral an triangular faces
        walln(SideN(s,n)) = 0.0
      end do
    end if
  end do 

  maxdis=0.0 
  do n=1,NN
    maxdis=max(walln(n),maxdis)
  end do

  do n=1,NN
    walln(n)=walln(n)/maxdis
  end do

  !------------------------------------------------------------!
  !   Calculate the interpolation factors for the cell sides   !
  !------------------------------------------------------------!
  if(rrun) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
  
      ! First cell
      xc1=xc(c1)
      yc1=yc(c1)
      zc1=zc(c1)
      dsc1=Distance(xc1,yc1,zc1,xsp(s), ysp(s), zsp(s))

      ! Second cell (pls. check if xsi=xc on the boundary)
      xc2=xc(c2)+Dx(s)
      yc2=yc(c2)+Dy(s)
      zc2=zc(c2)+Dz(s)
      dsc2=Distance(xc2,yc2,zc2,xsp(s), ysp(s), zsp(s))
  
      ! Interpolation factor
      f(s) = dsc2 / (dsc1+dsc2)   ! not checked
    end do 
  end if 

  end subroutine Compute_Grid_Geometry
