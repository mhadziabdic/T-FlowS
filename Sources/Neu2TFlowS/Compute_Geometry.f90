!==============================================================================!
  subroutine Compute_Geometry()
!------------------------------------------------------------------------------!
!   Calculates geometrical quantities of the grid.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
  real :: Distance_Squared    
  real :: Tol
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, n, s, ss, cc2, c_max, nnn, hh, mm
  integer              :: c11, c12, c21, c22, s1, s2, bou_cen
  integer              :: type_per, n_per, number_sides, dir, option
  integer              :: wall_mark, rot_dir, dir_face
  real                 :: xt(4), yt(4), zt(4), angle_face
  real                 :: xs2, ys2, zs2, x_a, y_a, z_a, x_b, y_b, z_b
  real                 :: x_c, y_c, z_c, Det
  real                 :: ab_i, ab_j, ab_k, ac_i, ac_j, ac_k, p_i, p_j, p_k
  real                 :: dsc1, dsc2, per_min, per_max
  real                 :: t, SurTot, angle 
  real                 :: xc1, yc1, zc1, xc2, yc2, zc2 
  real                 :: max_dis, tot_vol, min_vol, max_vol
  real                 :: xmin, xmax, ymin, ymax, zmin, zmax 
  real, allocatable    :: xnr(:), ynr(:), znr(:), xspr(:), yspr(:), zspr(:)
  real, allocatable    :: b_coor(:), phi_face(:)
  integer, allocatable :: b_face(:)
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
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
!     Sx * x_node + Sy * y_node + Sz * z_node = 
!     Sx * xsp(s) + Sy * ysp(s) + Sz * zsp(s)
!  
!     and the intersection is at:
!  
!         Sx*(xsp(s)-xc(c1)) + Sy*(ysp(s)-yc(c1) + Sz*(zsp(s)-zc(c1)) 
!     t = -----------------------------------------------------------
!                           rx*Sx + ry*Sy + rz*Sz
!  
!------------------------------------------------------------------------------!
 
  !-----------------------------------------!
  !   Calculate the cell centers            !
  !-----------------------------------------!
  !   => depends on: x_node,y_node,z_node   !
  !   <= gives:      xc,yc,zc c>0           !
  !-----------------------------------------!
  allocate(xc(-NbC:NC)); xc=0.0
  allocate(yc(-NbC:NC)); yc=0.0
  allocate(zc(-NbC:NC)); zc=0.0
  allocate(volume(NC)); volume=0.0

  do c=1,NC
    do n=1,CellN(c,0)
      xc(c) = xc(c) + x_node(CellN(c,n))/(1.0*CellN(c,0))
      yc(c) = yc(c) + y_node(CellN(c,n))/(1.0*CellN(c,0))
      zc(c) = zc(c) + z_node(CellN(c,n))/(1.0*CellN(c,0))
    end do
  end do

  write(*,*) '# Cell centers calculated !'

  !-----------------------------------------------------!
  !   Calculate:                                        ! 
  !      components of cell sides, cell side centers.   !
  !-----------------------------------------------------!
  !   => depends on: x_node,y_node,z_node               !
  !   <= gives:      Sx,Sy,Sz,xsp,yzp,zsp               !
  !-----------------------------------------------------!
  allocate(Sx(NS+max(NC,NBC))); Sx=0.0
  allocate(Sy(NS+max(NC,NBC))); Sy=0.0
  allocate(Sz(NS+max(NC,NBC))); Sz=0.0
  allocate(xsp(NS+max(NC,NBC))); xsp=0.0
  allocate(ysp(NS+max(NC,NBC))); ysp=0.0
  allocate(zsp(NS+max(NC,NBC))); zsp=0.0
  allocate(Dx(NS+max(NC,NBC))); Dx=0.0
  allocate(Dy(NS+max(NC,NBC))); Dy=0.0
  allocate(Dz(NS+max(NC,NBC))); Dz=0.0

  do s=1,NS
    do n=1,SideN(s,0)    ! for quadrilateral an triangular faces
      xt(n)=x_node(SideN(s,n))
      yt(n)=y_node(SideN(s,n))
      zt(n)=z_node(SideN(s,n))
    end do                       

    ! Cell side components
    if( SideN(s,0)  ==  4 ) then
      Sx(s)= 0.5 * ( (yt(2)-yt(1))*(zt(2)+zt(1))   &
                    +(yt(3)-yt(2))*(zt(2)+zt(3))   &
                    +(yt(4)-yt(3))*(zt(3)+zt(4))   &
                    +(yt(1)-yt(4))*(zt(4)+zt(1)) )
      Sy(s)= 0.5 * ( (zt(2)-zt(1))*(xt(2)+xt(1))   &
                    +(zt(3)-zt(2))*(xt(2)+xt(3))   &
                    +(zt(4)-zt(3))*(xt(3)+xt(4))   &
                    +(zt(1)-zt(4))*(xt(4)+xt(1)) )
      Sz(s)= 0.5 * ( (xt(2)-xt(1))*(yt(2)+yt(1))   & 
                    +(xt(3)-xt(2))*(yt(2)+yt(3))   &
                    +(xt(4)-xt(3))*(yt(3)+yt(4))   &
                    +(xt(1)-xt(4))*(yt(4)+yt(1)) )
    else if( SideN(s,0)  ==  3 ) then 
      Sx(s)= 0.5 * ( (yt(2)-yt(1))*(zt(2)+zt(1))   & 
                    +(yt(3)-yt(2))*(zt(2)+zt(3))   &
                    +(yt(1)-yt(3))*(zt(3)+zt(1)) )
      Sy(s)= 0.5 * ( (zt(2)-zt(1))*(xt(2)+xt(1))   &
                    +(zt(3)-zt(2))*(xt(2)+xt(3))   & 
                    +(zt(1)-zt(3))*(xt(3)+xt(1)) )
      Sz(s)= 0.5 * ( (xt(2)-xt(1))*(yt(2)+yt(1))   &
                    +(xt(3)-xt(2))*(yt(2)+yt(3))   & 
                    +(xt(1)-xt(3))*(yt(3)+yt(1)) )
    else
      write(*,*) 'calc4: something horrible has happened !'
      stop
    end if

    ! Barycenters
    if(SideN(s,0) == 4) then  
      xsp(s) = (xt(1)+xt(2)+xt(3)+xt(4))/4.0
      ysp(s) = (yt(1)+yt(2)+yt(3)+yt(4))/4.0
      zsp(s) = (zt(1)+zt(2)+zt(3)+zt(4))/4.0
    else if(SideN(s,0) == 3) then  
      xsp(s) = (xt(1)+xt(2)+xt(3))/3.0
      ysp(s) = (yt(1)+yt(2)+yt(3))/3.0
      zsp(s) = (zt(1)+zt(2)+zt(3))/3.0
    end if 

  end do ! through sides

  write(*,*) '# Cell face components calculated !'

  !--------------------------------------!
  !   Calculate boundary cell centers    !
  !--------------------------------------!
  !   => depends on: xc,yc,zc,Sx,Sy,Sz   !
  !   <= gives:      xc,yc,zc for c<0    !   
  !--------------------------------------!
  write(*,*) '#===================================='
  write(*,*) '# Position the boundary cell centres:'
  write(*,*) '#------------------------------------'
  write(*,*) '# Type 1 for barycentric placement'
  write(*,*) '# Type 2 for orthogonal placement'
  write(*,*) '#------------------------------------'
  read(*,*) bou_cen 

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    SurTot = sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))

    if(c2  < 0) then
      t = (   Sx(s)*(xsp(s)-xc(c1))                               &
            + Sy(s)*(ysp(s)-yc(c1))                               &
            + Sz(s)*(zsp(s)-zc(c1)) ) / SurTot
      xc(c2) = xc(c1) + Sx(s)*t / SurTot
      yc(c2) = yc(c1) + Sy(s)*t / SurTot
      zc(c2) = zc(c1) + Sz(s)*t / SurTot
      if(bou_cen == 1) then
        xc(c2) = xsp(s)
        yc(c2) = ysp(s)
        zc(c2) = zsp(s)
      end if
    endif 
  end do ! through sides

  !---------------------------------------------------------------!
  !   Move the centers of co-planar molecules towards the walls   !
  !---------------------------------------------------------------!
  !   => depends on: xc,yc,zc                                     !
  !   <= gives:      xc,yc,zc                                     !
  !   +  uses:       Dx,Dy,Dz                                     !
  !---------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2 > 0) then
      Dx(c1) = max( Dx(c1), abs( xc(c2)-xc(c1) ) )
      Dy(c1) = max( Dy(c1), abs( yc(c2)-yc(c1) ) )
      Dz(c1) = max( Dz(c1), abs( zc(c2)-zc(c1) ) )
      Dx(c2) = max( Dx(c2), abs( xc(c2)-xc(c1) ) )
      Dy(c2) = max( Dy(c2), abs( yc(c2)-yc(c1) ) )
      Dz(c2) = max( Dz(c2), abs( zc(c2)-zc(c1) ) )
    end if 
  end do ! through sides

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    if(c2 < 0) then
      if( Approx(Dx(c1), 0.0, 1.e-6) ) xc(c1) = 0.75*xc(c1) + 0.25*xc(c2) 
      if( Approx(Dy(c1), 0.0, 1.e-6) ) yc(c1) = 0.75*yc(c1) + 0.25*yc(c2) 
      if( Approx(Dz(c1), 0.0, 1.e-6) ) zc(c1) = 0.75*zc(c1) + 0.25*zc(c2) 
    end if 
  end do ! through sides

  Dx = 0.0
  Dy = 0.0
  Dz = 0.0

  !--------------------------------------------!
  !   Find the sides on the periodic boundary  !
  !--------------------------------------------!
  !   => depends on: xc,yc,zc,Sx,Sy,Sz         !
  !   <= gives:      Dx,Dy,Dz                  !
  !--------------------------------------------!
  allocate(b_coor(NS)); b_coor=0.0
  allocate(b_face(NS)); b_face=0

  !--------------------------------------------------------!
  !                                                        !
  !   Phase I  ->  find the sides on periodic boundaries   !
  !                                                        !
  !--------------------------------------------------------!
! write(*,*) 'Type 1 for fast but unreliable algorithm for periodic cells search.'
! write(*,*) 'Type 2 for slow but reliable algorithm for periodic cells search.'
! read(*,*) option

  option = 2

2 n_per = 0 
  write(*,*) '#======================================================'
  write(*,*) '# Type the periodic-boundary-condition marker number.'
  write(*,*) '# Type 0 if there is none !'
  write(*,*) '#------------------------------------------------------'
  write(*,*) '# (Please note that the periodic boundaries have to be' 
  write(*,*) '#  the last on the list of the boundary conditions.'
  write(*,*) '#  Their BC markers have to be larger than the markers'
  write(*,*) '#  of all the other boundary conditions.)'
  write(*,*) '#------------------------------------------------------'
  read(*,*) type_per
  if( type_per == 0 ) goto 1  

  if(option == 2) then

    write(*,*) '#========================================================'
    write(*,*) '# Insert the periodic direction (1 -> x, 2 -> y, 3 -> z)'
    write(*,*) '#--------------------------------------------------------'
    read(*,*) dir 
    write(*,*) 
    write(*,*) '#=============================================================='
    write(*,*) '# Enter the angle for the rotation of the coordiante system (in'
    write(*,*) '# degrees) and the axes of rotation (1 -> x, 2 -> y, 3 -> z)'
    write(*,*) '#--------------------------------------------------------------'
    write(*,*) '# (If the periodic direction is not parallel to the Caresian '
    write(*,*) '#  axis (x, y and z), the coordinate system has to be rotated'
    write(*,*) '#  in 2D'                                                     
    write(*,*) '#--------------------------------------------------------------'
    read(*,*) angle, rot_dir 

    angle = angle * PI / 180.0

    if(dir == 1) then
      x_a = 0.0
      y_a = 0.0
      z_a = 0.0       
      x_b = 0.0
      y_b = 1.0
      z_b = 0.0       
      x_c = 0.0
      y_c = 0.0
      z_c = 1.0       
    else if(dir == 2) then
      x_a = 0.0
      y_a = 0.0
      z_a = 0.0       
      x_b = 1.0
      y_b = 0.0
      z_b = 0.0       
      x_c = 0.0
      y_c = 0.0
      z_c = 1.0       
    else if(dir == 3) then
      x_a = 0.0
      y_a = 0.0
      z_a = 0.0       
      x_b = 1.0
      y_b = 0.0
      z_b = 0.0       
      x_c = 0.0
      y_c = 1.0
      z_c = 0.0       
    end if        

!   write(*,*) 
!   write(*,*) '# Enter the coordinates of three points that define'
!   write(*,*) '# periodic plane'
!   write(*,*) 
!   write(*,*) 'Point 1:'
!   read(*,*) x_a, y_a, z_a  !angle 
!   write(*,*) 
!   write(*,*) 'Point 2:'
!   read(*,*) x_b, y_b, z_b  !angle 
!   write(*,*) 
!   write(*,*) 'Point 3:'
!   read(*,*) x_c, y_c, z_c  !angle 

!   write(*,*) 
!   write(*,*) 'Enter approximative distance between the periodic faces'
!   read(*,*) per_max


    ab_i = x_b - x_a 
    ab_j = y_b - y_a 
    ab_k = z_b - z_a 

    ac_i = x_c - x_a 
    ac_j = y_c - y_a 
    ac_k = z_c - z_a

    p_i =  ab_j*ac_k - ac_j*ab_k 
    p_j = -ab_i*ac_k + ac_i*ab_k 
    p_k =  ab_i*ac_j - ac_i*ab_j

    angle_face = angle_face * PI / 180.0
 
    allocate(phi_face(NS)); phi_face=0.0
    allocate(xspr(NS)); xspr=0.0
    allocate(yspr(NS)); yspr=0.0
    allocate(zspr(NS)); zspr=0.0

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == type_per) then
          if( Approx(angle, 0.0, 1.e-6) ) then
            xspr(s) = xsp(s)  
            yspr(s) = ysp(s)  
            zspr(s) = zsp(s)  
          else 
            if(rot_dir==3) then
              xspr(s) = xsp(s)*cos(angle) + ysp(s)*sin(angle) 
              yspr(s) =-xsp(s)*sin(angle) + ysp(s)*cos(angle)
            else if(rot_dir==2) then
              xspr(s) = xsp(s)*cos(angle) + zsp(s)*sin(angle) 
              zspr(s) =-xsp(s)*sin(angle) + zsp(s)*cos(angle)
            else if(rot_dir==1) then
              yspr(s) = ysp(s)*cos(angle) + zsp(s)*sin(angle) 
              zspr(s) =-ysp(s)*sin(angle) + zsp(s)*cos(angle)
            end if
          end if
        end if
      end if
    end do
  end if

  xmin = +HUGE
  ymin = +HUGE
  zmin = +HUGE
  xmax = -HUGE
  ymax = -HUGE
  zmax = -HUGE


  b_coor=0.0
  b_face=0

  c = 0

  if(option == 1) then 
    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == type_per) then
          c = c + 1
          if(dir==1) b_coor(c) = xsp(s)*1000000.0 + ysp(s)*10000.0 + zsp(s)
          if(dir==2) b_coor(c) = ysp(s)*1000000.0 + xsp(s)*10000.0 + zsp(s)
          if(dir==3) b_coor(c) = zsp(s)*1000000.0 + xsp(s)*10000.0 + ysp(s)
          b_face(c) = s
        end if
      end if
    end do
    call Sort_Double_Carry_Int(b_coor,b_face,c,2)
  else if(option == 2) then 
    c_max = 0
    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == type_per) then
          c_max = c_max + 1
        end if 
      end if 
    end do
    Tol = 0.0001 
    
10 continue


    nnn = 0
    hh = 0 
    mm = 0 
    c = 0
  
    per_max = -HUGE
    per_min =  HUGE 

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == type_per) then
          Det = (p_i*(xsp(s)) + p_j*(ysp(s)) + p_k*(zsp(s)))/sqrt(p_i*p_i + p_j*p_j + p_k*p_k)
          per_min = min(per_min, Det)          
          per_max = max(per_max, Det)          
        end if
      end if
    end do
    
    per_max = 0.5*(per_max + per_min)

    do s=1,NS
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == type_per) then
          c = c + 1
          Det = (p_i*(xsp(s)) + p_j*(ysp(s)) + p_k*(zsp(s)))  &
              / sqrt(p_i*p_i + p_j*p_j + p_k*p_k)

          if(dir==1) then
            if((Det) < (per_max)) then
              hh = hh + 1
              b_coor(hh) = hh
              b_face(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == type_per) then 
                    Det = (p_i*(xsp(ss)) + p_j*(ysp(ss)) + p_k*(zsp(ss)))  &
                        / sqrt(p_i*p_i + p_j*p_j + p_k*p_k)
                    if((Det) > (per_max)) then
                      if((abs(zsp(ss)-zsp(s))) < Tol.and.         &
                         (abs(yspr(ss)-yspr(s))) < Tol) then
                         mm = hh + c_max/2 
                         b_coor(mm) = mm
                         b_face(mm) = ss
                         nnn = nnn + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if


          if(dir==2) then
            if((Det) < (per_max)) then
              hh = hh + 1
              b_coor(hh) = hh
              b_face(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == type_per) then 
                    Det = (p_i*(xsp(ss)) + p_j*(ysp(ss)) + p_k*(zsp(ss))) &
                        / sqrt(p_i*p_i + p_j*p_j + p_k*p_k)
                    if((Det) > (per_max)) then
                      if(abs((zsp(ss)-zsp(s))) < Tol.and.         &
                        abs((xspr(ss)-xspr(s))) < Tol) then
                        mm = hh + c_max/2 
                        b_coor(mm) = mm
                        b_face(mm) = ss
                        nnn = nnn + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if

          if(dir==3) then
            if((Det) < (per_max)) then
              hh = hh + 1
              b_coor(hh) = hh
              b_face(hh) = s
              do ss=1,NS
                cc2 = SideC(2,ss)
                if(cc2 < 0) then
                  if(BCmark(cc2) == type_per) then
                    Det = (p_i*(xsp(ss)) + p_j*(ysp(ss)) + p_k*(zsp(ss)))  &
                        / sqrt(p_i*p_i + p_j*p_j + p_k*p_k)
                    if((Det) > (per_max)) then
                      if(abs((xsp(ss)-xsp(s))) < Tol.and.         &
                        abs((ysp(ss)-ysp(s))) < Tol) then
                        mm = hh + c_max/2 
                        b_coor(mm) = mm
                        b_face(mm) = ss
                        nnn = nnn + 1
                      end if 
                    end if 
                  end if
                end if  
              end do
            end if          
          end if

        end if 
      end if 
    end do

    write(*,*)'Iterating search for periodic cells: ',  &
    'Target: ', c_max/2, 'Result: ',nnn, 'Tolerance: ',Tol

    if(nnn == c_max/2) then
      continue
    else
      Tol = Tol*0.5 
      go to 10
    end if

    deallocate(phi_face)
    deallocate(xspr)
    deallocate(yspr)
    deallocate(zspr)

    call Sort_Double_Carry_Int(b_coor,b_face,c,2)
  end if  ! end option

  do s=1,c/2
    s1 = b_face(s)
    s2 = b_face(s+c/2)
    c11 = SideC(1,s1)  ! cell 1 for side 1
    c21 = SideC(2,s1)  ! cell 2 for cell 1
    c12 = SideC(1,s2)  ! cell 1 for side 2
    c22 = SideC(2,s2)  ! cell 2 for side 2
    SideC(0,s1) = s2   ! just to remember where it was coppied from
    SideC(2,s1) = c12 
    SideC(1,s2) = 0    ! c21 
    SideC(2,s2) = 0    ! c21 
  end do

  n_per = c/2
  write(*,*) 'Phase I: periodic cells: ', n_per
  go to 2

  !----------------------------------------------------!
  !                                                    !
  !   Phase II  ->  similar to the loop in Generator   !
  !                                                    !
  !----------------------------------------------------!
1 continue 
  NSsh = 0
  do s=1,NS

    ! Initialize
    Dx(s)=0.0
    Dy(s)=0.0
    Dz(s)=0.0

    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2   >  0) then

      ! Scalar product of the side with line c1-c2 is good criteria
      if( (Sx(s) * (xc(c2)-xc(c1) )+                              &
           Sy(s) * (yc(c2)-yc(c1) )+                              &
           Sz(s) * (zc(c2)-zc(c1) ))  < 0.0 ) then

        NSsh = NSsh + 2

        ! Find the coordinates of ...
        if(SideN(s,0) == 4) then

          ! Coordinates of the shadow face
          xs2=xsp(SideC(0,s))
          ys2=ysp(SideC(0,s))
          zs2=zsp(SideC(0,s))

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
          SideN(NS+NSsh,1) = SideN(SideC(0,s),1) 
          SideN(NS+NSsh,2) = SideN(SideC(0,s),2)
          SideN(NS+NSsh,3) = SideN(SideC(0,s),3)
          SideN(NS+NSsh,4) = SideN(SideC(0,s),4)
          Sx(NS+NSsh) = Sx(s)
          Sy(NS+NSsh) = Sy(s)
          Sz(NS+NSsh) = Sz(s)
          xsp(NS+NSsh) = xs2
          ysp(NS+NSsh) = ys2
          zsp(NS+NSsh) = zs2
        else if(SideN(s,0) == 3) then

          ! Coordinates of the shadow face
          xs2=xsp(SideC(0,s))
          ys2=ysp(SideC(0,s))
          zs2=zsp(SideC(0,s))
 
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
          SideN(NS+NSsh,1) = SideN(SideC(0,s),1) 
          SideN(NS+NSsh,2) = SideN(SideC(0,s),2)
          SideN(NS+NSsh,3) = SideN(SideC(0,s),3)
          Sx(NS+NSsh) = Sx(s)
          Sy(NS+NSsh) = Sy(s)
          Sz(NS+NSsh) = Sz(s)
          xsp(NS+NSsh) = xs2
          ysp(NS+NSsh) = ys2
          zsp(NS+NSsh) = zs2
        end if

        Dx(s)=xsp(s)-xs2  !
        Dy(s)=ysp(s)-ys2  ! later: xc2 = xc2 + Dx  
        Dz(s)=zsp(s)-zs2  !

      endif !  S*(c2-c1) < 0.0
    end if  !  c2 > 0
  end do    !  sides
  write(*,*) 'Phase II: number of shadow faces: ', NSsh

  !-------------------------------------------------------!
  !                                                       !
  !   Phase III  ->  find the new numbers of cell faces   !
  !                                                       !
  !-------------------------------------------------------!
  number_sides = 0
  do s=1,NS+NSsh
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if(c1 > 0) then
      number_sides = number_sides  + 1
      NewS(s) = number_sides 
    else
      NewS(s) = -1
    end if
  end do
  write(*,'(A21,I9,Z9)') '# Old number of sides: ', NS, NS
  write(*,'(A21,I9,Z9)') '# New number of sides: ', &
                          number_sides-NSsh,number_sides-NSsh
  
  !--------------------------------------!
  !                                      !
  !   Phase IV  ->  compress the sides   !
  !                                      !
  !--------------------------------------!
  do s=1,NS+NSsh
    if(NewS(s) > 0) then
      SideC(1,NewS(s)) = SideC(1,s) 
      SideC(2,NewS(s)) = SideC(2,s)
      SideN(NewS(s),0) = SideN(s,0)
      SideN(NewS(s),1) = SideN(s,1)
      SideN(NewS(s),2) = SideN(s,2)
      SideN(NewS(s),3) = SideN(s,3)
      SideN(NewS(s),4) = SideN(s,4)
      xsp(NewS(s)) = xsp(s)
      ysp(NewS(s)) = ysp(s)
      zsp(NewS(s)) = zsp(s)
      Sx(NewS(s)) = Sx(s)
      Sy(NewS(s)) = Sy(s)
      Sz(NewS(s)) = Sz(s)
      Dx(NewS(s)) = Dx(s)
      Dy(NewS(s)) = Dy(s)
      Dz(NewS(s)) = Dz(s)
    end if
  end do 
  NS = number_sides-NSsh

  !-----------------------------------!
  !   Check the periodic boundaries   !
  !-----------------------------------!
  max_dis = 0.0 
  do s=1,NS-NSsh
    max_dis = max(max_dis, (Dx(s)*Dx(s)+Dy(s)*Dy(s)+Dz(s)*Dz(s)))
  end do
  write(*,*) '# Maximal distance of periodic boundary is:', sqrt(max_dis)

  !----------------------------------!
  !   Calculate the cell volumes     !
  !----------------------------------!
  !   => depends on: xc,yc,zc,       !
  !                  Dx,Dy,Dz,       !
  !                  xsp, ysp, zsp   !
  !   <= gives:      volume          !
  !----------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)   

    volume(c1) = volume(c1) + xsp(s)*Sx(s)  &
                            + ysp(s)*Sy(s)  &
                            + zsp(s)*Sz(s)
    if(c2 > 0) then
      volume(c2) = volume(c2) - (xsp(s)-Dx(s))*Sx(s)  &
                              - (ysp(s)-Dy(s))*Sy(s)  &
                              - (zsp(s)-Dz(s))*Sz(s)
    end if
  end do
  volume = volume/3.0
  c1 = 0
  min_vol =  1E+30
  max_vol = -1E+30
  tot_vol = 0.0
  do c=1,NC
    tot_vol = tot_vol + volume(c)
    min_vol = min(min_vol, volume(c))
    max_vol = max(max_vol, volume(c))
  end do
  write(*,*) '# Minimal cell volume is: ', min_vol
  write(*,*) '# Maximal cell volume is: ', max_vol
  write(*,*) '# Total domain volume is: ', tot_vol
  write(*,*) '# Cell volumes calculated !'

  if(min_vol < 0.0) then
    write(*,*) '# Negative volume occured! Another, slower, algoritham should be run !'
    write(*,*) '# Execution will be halt now! '
    stop
  end if 
 
  deallocate(b_coor)
  deallocate(b_face)
 

    !------------------------------------------!
    !     Calculate delta                      !
    !------------------------------------------!
    !     => depends on: x_node,y_node,z_node  !
    !     <= gives:      delta                 !
    !------------------------------------------!
    allocate(delta(-NbC:NC));  delta=0.0

    do c=1,NC
      delta(c)=0.0
      xmin = +HUGE
      ymin = +HUGE
      zmin = +HUGE
      xmax = -HUGE
      ymax = -HUGE
      zmax = -HUGE
      do n=1,CellN(c,0)
        xmin = min(xmin, x_node(CellN(c,n)))
        ymin = min(ymin, y_node(CellN(c,n)))
        zmin = min(zmin, z_node(CellN(c,n)))
        xmax = max(xmax, x_node(CellN(c,n)))
        ymax = max(ymax, y_node(CellN(c,n)))
        zmax = max(zmax, z_node(CellN(c,n)))
      end do
      delta(c) = xmax-xmin
      delta(c) = max(delta(c), (ymax-ymin))
      delta(c) = max(delta(c), (zmax-zmin))
    end do

  !------------------------------------------------------------------!
  !   Calculate distance from the cell center to the nearest wall.   !
  !------------------------------------------------------------------!
  !     => depends on: xc,yc,zc inside and on the boundary.          !
  !     <= gives:      WallDs i                                      !
  !------------------------------------------------------------------!
  allocate(WallDs(-NbC:NC)); WallDs = HUGE

  write(*,*) '#================================================================'
  write(*,*) '# Type the total number of wall boundary conditions:'
  write(*,*) '#----------------------------------------------------------------'
  write(*,*) '# (Please note that the walls have to be the first on the list)'
  write(*,*) '# of the boundary conditions. Their BC markers have to be smaller'
  write(*,*) '# than the markers of the other boundary conditions.)'
  write(*,*) '#----------------------------------------------------------------'
  read(*,*) wall_mark
 
  if(wall_mark == 0) then
    WallDs = 1.0
    write(*,*) '# Distance to the wall set to 1.0 everywhere !'
  else
    do c1=1,NC
      if(mod(c1,10000) == 0) then
        write(*,'(a2, f5.0, a14)') ' #', (100.*c1/(1.*NC)), ' % complete...'
      endif
      do c2=-1,-NbC,-1
        if(BCmark(c2) <= wall_mark) then
          WallDs(c1)=min(WallDs(c1),                       &
          Distance_Squared(xc(c1),yc(c1),zc(c1),xc(c2),yc(c2),zc(c2)))
        end if
      end do
    end do

    WallDs = sqrt(WallDs)

    write(*,*) '# Distance to the wall calculated !'
  end if


  !------------------------------------------------------------!
  !   Calculate the interpolation factors for the cell sides   !
  !------------------------------------------------------------!
  allocate(f(NS+max(NC,NBC))); f=0.0          

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!----- first cell
    xc1  = xc(c1)
    yc1  = yc(c1)
    zc1  = zc(c1)
    dsc1 = Distance(xc1, yc1, zc1, xsp(s), ysp(s), zsp(s))

!----- second cell (pls. check if xsi=xc on the boundary)
    xc2  = xc(c2)+Dx(s)
    yc2  = yc(c2)+Dy(s)
    zc2  = zc(c2)+Dz(s)
    dsc2 = Distance(xc2, yc2, zc2, xsp(s), ysp(s), zsp(s))

!----- interpolation factor
    f(s) = dsc2 / (dsc1+dsc2)   ! not checked
  end do 

  write(*,*) '# Interpolation factors calculated !'


  return

3 write(*,*) '# Horror ! Negative volume between cells ', c1, ' and ', c2
  stop

  end subroutine Compute_Geometry
