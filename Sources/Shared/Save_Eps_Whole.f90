!==============================================================================!
  subroutine Save_Eps_Whole(NSsh0)
!------------------------------------------------------------------------------!
!   Saves the whole grid in encapsulated postscript.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!---------------------------------[Arguments]----------------------------------!
  integer :: NSsh0  ! if 0, shadow will not be drawn 
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n, s, s0, c1, c2
  integer             :: xmaxb, xminb, ymaxb, yminb, xlegend
  character(len=80)   :: name_eps, answer, colour
  real                :: sclf, sclp, xmax,xmin,ymax,ymin,zmax,zmin
  real                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,    &
                         xk,yk,zk,alfa,beta,gama,nx,ny,nz,shade, &
                         xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4 
  integer,allocatable :: indx(:)
  real,allocatable    :: work(:)
  real                :: red(10), green(10), blue(10)
  integer             :: ix1, ix2, iy1, iy2, boxsize, BCcount(0:11)
!==============================================================================!

  ! Allocate the memory
  allocate(indx(max(NS+NSsh0,MAXS))); indx=0
  allocate(work(max(NS+NSsh0,MAXS))); work=0

  BCcount = 0

  write(*,*) '#============================================='
  write(*,*) '# Making a 3D shaded .eps figure of the domain'
  write(*,*) '#---------------------------------------------'

  !---------------------------------------!
  !   Set the boundary condition colors   !
  !---------------------------------------!
  n=1  ! colour count

  red(n) = 1.00; green(n) = 0.00; blue(n) = 0.00;  n=n+1 ! red
  red(n) = 0.00; green(n) = 1.00; blue(n) = 0.00;  n=n+1 ! green
  red(n) = 0.00; green(n) = 0.00; blue(n) = 1.00;  n=n+1 ! blue
  red(n) = 1.00; green(n) = 0.00; blue(n) = 1.00;  n=n+1 ! magenta
  red(n) = 1.00; green(n) = 1.00; blue(n) = 0.00;  n=n+1 ! yellow
  red(n) = 0.00; green(n) = 1.00; blue(n) = 1.00;  n=n+1 ! cyan
  red(n) = 1.00; green(n) = 0.50; blue(n) = 0.50;  n=n+1 ! light Red
  red(n) = 0.50; green(n) = 1.00; blue(n) = 0.50;  n=n+1 ! light Green
  red(n) = 0.50; green(n) = 0.50; blue(n) = 1.00;  n=n+1 ! light Blue
  red(n) = 0.50; green(n) = 0.50; blue(n) = 0.50;  n=n+1 ! gray           

  !------------------------------!
  !   Input camera coordinates   !
  !------------------------------!
1 write(6,*) '# Enter the camera coordinates (skip to exit): '
  write(*,*) '#---------------------------------------------'
  call ReadC(5,inp,tn,ts,te)
  if(tn == 1) then 
    read(inp, *) answer
    call To_Upper_Case(answer)
    if(answer == 'SKIP') return  
  else if(tn == 3) then 
    read(inp, *) xk, yk, zk  
  end if  
  alfa = acos( xk / sqrt(xk*xk+yk*yk) )
  beta = acos( yk / sqrt(xk*xk+yk*yk) )
  gama = acos( zk / sqrt(xk*xk+yk*yk+zk*zk) )

  nx = xk/sqrt(xk*xk+yk*yk+zk*zk)
  ny = yk/sqrt(xk*xk+yk*yk+zk*zk)
  nz = zk/sqrt(xk*xk+yk*yk+zk*zk)

  !----------------------!
  !                      !
  !   Create .eps file   !
  !                      !
  !----------------------!
  write(6,*) '# G-> Gray or C-> Coloured (by boundary conditions): '
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) colour 
  call To_Upper_Case(colour);
  write(6,*) '# Enter the file name (without extension): '
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) name_eps 
  name_eps(len_trim(name_eps)+1:len_trim(name_eps)+4) = '.eps'
  write(6, *) '# Now creating the file:', name_eps

  xmax=maxval(x_node(1:NN))
  ymax=maxval(y_node(1:NN))
  zmax=maxval(z_node(1:NN))
  xmin=minval(x_node(1:NN))
  ymin=minval(y_node(1:NN))
  zmin=minval(z_node(1:NN))
  sclf = 100000.0/max((xmax-xmin),(ymax-ymin),(zmax-zmin))
  sclp = 0.005 

!---- Find the bounding box
  xminb=+1000000
  xmaxb=-1000000
  yminb=+1000000
  ymaxb=-1000000
  do n=1,Nn
    if(xk  < 0.0 .and. yk  > 0.0) then
      xp1=-x_node(n)*sin(alfa)-y_node(n)*sin(beta)
    else if(xk  > 0.0 .and. yk  < 0.0) then
      xp1=x_node(n)*sin(alfa)+y_node(n)*sin(beta)
    else if(xk  > 0.0 .and. yk  > 0.0) then
      xp1=-x_node(n)*sin(alfa)+y_node(n)*sin(beta)
    else
      xp1=x_node(n)*sin(alfa)-y_node(n)*sin(beta)
    end if
    xp1=xp1*sclf*sclp
    xmaxb=max(xmaxb,int(xp1))
    xminb=min(xminb,int(xp1))

    yp1=(-x_node(n)*cos(alfa)-y_node(n)*cos(beta))*cos(gama)+z_node(n)*sin(gama) 
    yp1=yp1*sclf*sclp
    ymaxb=max(ymaxb,int(yp1))
    yminb=min(yminb,int(yp1))
  end do
  boxsize = (ymaxb-yminb) / 20
  if(colour(1:1) == 'C') then
    xlegend = boxsize * 8  
  else
    xlegend = 0
  end if

  open(9, file=name_eps)

  write(*,*) 'File opened'

  write(9, '(A24)') '%!PS-Adobe-2.0 EPSF-1.2 '
  write(9, '(A24)') '%% Created by:  TFlowS %%'
  write(9, '(A15,4I7)') '%%BoundingBox: ',                          &
                        xminb-2,yminb-2,xmaxb+2+xlegend,ymaxb+2 
  write(9, '(A24)') '% Abrevations:          '
  write(9, '(A24)') '/cp   {closepath}    def'
  write(9, '(A24)') '/f    {fill}         def'
  write(9, '(A24)') '/gr   {grestore}     def'
  write(9, '(A24)') '/gs   {gsave}        def'
  write(9, '(A24)') '/l    {lineto}       def'
  write(9, '(A24)') '/m    {moveto}       def'
  write(9, '(A24)') '/np   {newpath}      def'
  write(9, '(A24)') '/s    {stroke}       def'
  write(9, '(A24)') '/sc   {scale}        def'
  write(9, '(A24)') '/sg   {setgray}      def'
  write(9, '(A24)') '/sh   {show}         def'
  write(9, '(A24)') '/slw  {setlinewidth} def'
  write(9, '(A24)') '/srgb {setrgbcolor}  def'
  write(9, '(A24)') '% End of brevations.    '

  write(9,'(A4)') ' gs '
  write(9, '(A24)') '% Scale:                '
  write(9, '(2F10.6,A4)') sclp,sclp, ' sc '

  write(9, '(A)') '%% Font'
  write(9, '(A)') '/Arial findfont'
  write(9, '(I6,A)') int(real(boxsize)/sclp), ' scalefont'
  write(9, '(A)') 'setfont'

  do s=1,NS+NSsh0
    indx(s) = s
    work(s)  = Distance(xk,yk,zk,xsp(s),ysp(s),zsp(s))
  end do
  call Sort_Real_By_Index(work,indx,NS+NSsh0,-2)

  do s0=1,NS+NSsh0
    s=indx(s0)

    shade=                                                          &
   (Sx(s)*nx+Sy(s)*ny+Sz(s)*nz)                                     &
    /sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))
    shade=abs(shade)
    shade=0.4+0.6*shade

    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if(c2 < 0 .or. material(c1) /= material(c2) ) then 

      BCcount(BCmark(c2)) = BCcount(BCmark(c2)) + 1

      x1 = x_node(SideN(s,1))
      x2 = x_node(SideN(s,2))
      x3 = x_node(SideN(s,3))
      x4 = x_node(SideN(s,4))

      y1 = y_node(SideN(s,1))
      y2 = y_node(SideN(s,2))
      y3 = y_node(SideN(s,3))
      y4 = y_node(SideN(s,4))

      z1 = z_node(SideN(s,1))
      z2 = z_node(SideN(s,2))
      z3 = z_node(SideN(s,3))
      z4 = z_node(SideN(s,4))

      if(xk  < 0.0 .and. yk  > 0.0) then
        xp1=-x1*sin(alfa)-y1*sin(beta)
        xp2=-x2*sin(alfa)-y2*sin(beta)
        xp3=-x3*sin(alfa)-y3*sin(beta)
        xp4=-x4*sin(alfa)-y4*sin(beta)
      else if(xk  > 0.0 .and. yk  < 0.0) then
        xp1=x1*sin(alfa)+y1*sin(beta)
        xp2=x2*sin(alfa)+y2*sin(beta)
        xp3=x3*sin(alfa)+y3*sin(beta)
        xp4=x4*sin(alfa)+y4*sin(beta)
      else if(xk  > 0.0 .and. yk  > 0.0) then 
        xp1=-x1*sin(alfa)+y1*sin(beta)
        xp2=-x2*sin(alfa)+y2*sin(beta)
        xp3=-x3*sin(alfa)+y3*sin(beta)
        xp4=-x4*sin(alfa)+y4*sin(beta)
      else
        xp1=x1*sin(alfa)-y1*sin(beta)
        xp2=x2*sin(alfa)-y2*sin(beta)
        xp3=x3*sin(alfa)-y3*sin(beta)
        xp4=x4*sin(alfa)-y4*sin(beta)
      end if

      yp1=(-x1*cos(alfa)-y1*cos(beta))*cos(gama) + z1*sin(gama)
      yp2=(-x2*cos(alfa)-y2*cos(beta))*cos(gama) + z2*sin(gama)
      yp3=(-x3*cos(alfa)-y3*cos(beta))*cos(gama) + z3*sin(gama)
      yp4=(-x4*cos(alfa)-y4*cos(beta))*cos(gama) + z4*sin(gama)

      if(s <= NS) then
        if(colour(1:1) == 'G') then
          if(SideN(s,0) == 4)                                       &
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A9,F4.2,A15)')  &
                       'gs np 0 slw ',                                 &
                    int(sclf*xp1),int(sclf*yp1), ' m ',             &
                    int(sclf*xp2),int(sclf*yp2), ' l ',             &
                    int(sclf*xp3),int(sclf*yp3), ' l ',             &
                    int(sclf*xp4),int(sclf*yp4), ' l ',             &
                    ' cp   gs ',shade,' sg f gr   s gr'
          if(SideN(s,0) == 3)                                       &
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,A9,F4.2,A15)')         &
                      'gs np 0 slw ',                                 &
                    int(sclf*xp1),int(sclf*yp1), ' m ',             &
                    int(sclf*xp2),int(sclf*yp2), ' l ',             &
                    int(sclf*xp3),int(sclf*yp3), ' l ',             &
                    ' cp   gs ',shade,' sg f gr   s gr'
        else
          if(SideN(s,0) == 4)                                       &
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A7,3F5.2,A15)') &
                      'gs np 0 slw ',                                 &
                    int(sclf*xp1),int(sclf*yp1), ' m ',             & 
                    int(sclf*xp2),int(sclf*yp2), ' l ',             &
                    int(sclf*xp3),int(sclf*yp3), ' l ',             &
                    int(sclf*xp4),int(sclf*yp4), ' l ',             &
                    ' cp gs ',                                      &
                    red(BCmark(c2)),                                &
                    green(BCmark(c2)),                              &
                    blue(BCmark(c2)),                               &
                    ' srgb f gr s gr'

          if(SideN(s,0) == 3)                                       &
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,A7,3F5.2,A15)')        &
                    'gs np 0 slw ',                                 &
                    int(sclf*xp1),int(sclf*yp1), ' m ',             & 
                    int(sclf*xp2),int(sclf*yp2), ' l ',             &
                    int(sclf*xp3),int(sclf*yp3), ' l ',             &
                    ' cp gs ',                                      &
                    red(BCmark(c2)),                                &
                    green(BCmark(c2)),                              &
                    blue(BCmark(c2)),                               &
                    ' srgb f gr s gr'
        end if ! colour
      else if(s > NS) then
        if(SideN(s,0) == 4) then
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A9)')           &
                        'gs np 0 slw ',                                 &
                    int(sclf*xp1),int(sclf*yp1), ' m ',             &
                    int(sclf*xp2),int(sclf*yp2), ' l ',             &
                    int(sclf*xp3),int(sclf*yp3), ' l ',             &
                    int(sclf*xp4),int(sclf*yp4), ' l ',             &
                    ' cp s gr '
        else if(SideN(s,0) == 3) then
          write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,A9)')           &
                      'gs np 0 slw ',                          &
                    int(sclf*xp1),int(sclf*yp1), ' m ',      &
                    int(sclf*xp2),int(sclf*yp2), ' l ',      &
                    int(sclf*xp3),int(sclf*yp3), ' l ',      &
                    ' cp s gr '
        end if
      end if
    end if

  end do 

  !-----------------!
  !   Plot legend   !
  !-----------------!
  if(colour(1:1) == 'C') then
    ix1 = xmaxb + boxsize/2
    iy1 = ymaxb
    ix2 = ix1+boxsize
    iy2 = iy1+boxsize
    do n=1,19
      if( BCcount(n) > 0 ) then
        iy1 = iy1 - boxsize
        iy2 = iy2 - boxsize
        write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A7,3F5.2,A15)') &
                  'gs np 0 slw ',                                 &
                  int(1.0/sclp*ix1),int(1.0/sclp*iy1), ' m ',     &
                  int(1.0/sclp*ix1),int(1.0/sclp*iy2), ' l ',     &
                  int(1.0/sclp*ix2),int(1.0/sclp*iy2), ' l ',     &
                  int(1.0/sclp*ix2),int(1.0/sclp*iy1), ' l ',     &
                  ' cp gs ',                                      &
                  red(n), green(n), blue(n),                      &
                  ' srgb f gr s gr'
        write(9,'(A6,2I8,A13,I3, A1, A6)')                         &
                  'gs np ',                                        &
                  int(1.0/sclp*(ix2+boxsize/2)),int(1.0/sclp*iy1), &
                  ' m (boundary:', n, ')',                         &
                  ' sh gr'              
      end if
    end do
  end if

  write(9,'(A4)') ' gr '

  close(9)

  goto 1

  ! Free the memory
  deallocate(indx)
  deallocate(work)

  end subroutine Save_Eps_Whole
