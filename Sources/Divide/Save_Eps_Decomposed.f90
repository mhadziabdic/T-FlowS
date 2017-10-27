!==============================================================================!
  subroutine Save_Eps_Decomposed()
!------------------------------------------------------------------------------!
!   Saves the whole grid in encapsulated postscript.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n, s, s0, c1, c2
  integer             :: xmaxb, xminb, ymaxb, yminb, xlegend
  character(len=80)   :: name_eps, answer
  real                :: sclf, sclp, xmax,xmin,ymax,ymin,zmax,zmin
  real                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,        &
                         xk,yk,zk,alfa,beta,gama,nx,ny,nz,shade,     &
                         xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4 
  real                :: red(MAXPRO), green(MAXPRO), blue(MAXPRO)
  integer,allocatable :: indx(:)
  real,allocatable    :: work(:) 
  integer             :: ix1, ix2, iy1, iy2, boxsize
!==============================================================================!

  ! Allocate the memory
  allocate(indx(NS)); indx=0
  allocate(work(NS)); work=0   

  !------------------------------!
  !   Set the processor colors   !
  !------------------------------! 
  n=1  ! colour count

  ! Few basic colors
  red(n) = 1.00; green(n) = 0.00; blue(n) = 0.00;  n=n+1 ! red
  red(n) = 1.00; green(n) = 1.00; blue(n) = 0.00;  n=n+1 ! yellow
  red(n) = 0.00; green(n) = 1.00; blue(n) = 0.00;  n=n+1 ! green
  red(n) = 0.00; green(n) = 1.00; blue(n) = 1.00;  n=n+1 ! cyan
  red(n) = 0.00; green(n) = 0.00; blue(n) = 1.00;  n=n+1 ! blue
  red(n) = 1.00; green(n) = 0.00; blue(n) = 1.00;  n=n+1 ! magenta

  ! Few basic colors but lighter 
  red(n) = 1.00; green(n) = 0.50; blue(n) = 0.50;  n=n+1 ! light red
  red(n) = 1.00; green(n) = 1.00; blue(n) = 0.50;  n=n+1 ! light yellow
  red(n) = 0.50; green(n) = 1.00; blue(n) = 0.50;  n=n+1 ! light green
  red(n) = 0.50; green(n) = 1.00; blue(n) = 1.00;  n=n+1 ! light cyan
  red(n) = 0.50; green(n) = 0.50; blue(n) = 1.00;  n=n+1 ! light blue
  red(n) = 1.00; green(n) = 0.50; blue(n) = 1.00;  n=n+1 ! light magenta

  ! Few basic colors but darker
  red(n) = 0.50; green(n) = 0.00; blue(n) = 0.00;  n=n+1 ! dark red
  red(n) = 0.50; green(n) = 0.50; blue(n) = 0.00;  n=n+1 ! dark yellow
  red(n) = 0.00; green(n) = 0.50; blue(n) = 0.00;  n=n+1 ! dark green
  red(n) = 0.00; green(n) = 0.50; blue(n) = 0.50;  n=n+1 ! dark cyan
  red(n) = 0.00; green(n) = 0.00; blue(n) = 0.50;  n=n+1 ! dark blue
  red(n) = 0.50; green(n) = 0.00; blue(n) = 0.50;  n=n+1 ! dark magenta

  !------------------------------!
  !   Input camera coordinates   !
  !------------------------------!
1 write(6,*) 'Enter the camera coordinates (skip to exit): '
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
  write(6,*) 'Enter the file name (without extension): '
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) name_eps 
  name_eps(len_trim(name_eps)+1:len_trim(name_eps)+4) = '.eps'
  write(6, *) 'Now creating the file:', name_eps

  xmax=maxval(grid % nodes(1:NN) % x)
  ymax=maxval(grid % nodes(1:NN) % y)
  zmax=maxval(grid % nodes(1:NN) % z)
  xmin=minval(grid % nodes(1:NN) % x)
  ymin=minval(grid % nodes(1:NN) % y)
  zmin=minval(grid % nodes(1:NN) % z)
  sclf = 100000.0/max((xmax-xmin),(ymax-ymin),(zmax-zmin))
  sclp = 0.005 

  ! Find the bounding box
  xminb=+1000000
  xmaxb=-1000000
  yminb=+1000000
  ymaxb=-1000000
  do n=1,Nn
    if(xk  < 0.0 .and. yk  > 0.0) then
      xp1= - grid % nodes(n) % x * sin(alfa) - grid % nodes(n) % y *sin(beta)
    else if(xk  > 0.0 .and. yk  < 0.0) then
      xp1=   grid % nodes(n) % x * sin(alfa) + grid % nodes(n) % y *sin(beta)
    else if(xk  > 0.0 .and. yk  > 0.0) then
      xp1= - grid % nodes(n) % x * sin(alfa) + grid % nodes(n) % y *sin(beta)
    else
      xp1=   grid % nodes(n) % x * sin(alfa) - grid % nodes(n) % y *sin(beta)
    end if
    xp1=xp1*sclf*sclp
    xmaxb=max(xmaxb,int(xp1))
    xminb=min(xminb,int(xp1))

    yp1=( - grid % nodes(n) % x * cos(alfa)                  &
          - grid % nodes(n) % y * cos(beta) ) * cos(gama) +  &
            grid % nodes(n) % z * sin(gama) 
    yp1=yp1*sclf*sclp
    ymaxb=max(ymaxb,int(yp1))
    yminb=min(yminb,int(yp1))
  end do
  boxsize = (ymaxb - yminb) / 20
  xlegend = boxsize * 8 

  open(9, file=name_eps)

  ! Create the header of the .eps file
  write(9, '(A24)') '%!PS-Adobe-2.0 EPSF-1.2 '
  write(9, '(A24)') '%% Created by:  T-Rex %%'
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
  write(9, '(A24)') '% End of abrevations.   '

  write(9,'(A4)') ' gs '
  write(9, '(A24)') '% Scale:                '
  write(9, '(2F10.6,A4)') sclp,sclp, ' sc '

  write(9, '(A)') '%% Font'
  write(9, '(A)') '/Arial findfont'
  write(9, '(I6,A)') int(real(boxsize)/sclp), ' scalefont' 
  write(9, '(A)') 'setfont'

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    indx(s) = s
    work(s) = Distance(xk,yk,zk,                              &
                       f(s)*xc(c1)+(1.-f(s))*(xc(c2)+Dx(s)),  &
                       f(s)*yc(c1)+(1.-f(s))*(yc(c2)+Dy(s)),  &
                       f(s)*zc(c1)+(1.-f(s))*(zc(c2)+Dz(s)) )
  end do
  call Sort_Real_By_Index(work,indx,NS,-2)

  do s0=1,NS
    s=indx(s0)

    shade=                                                          &
   (Sx(s)*nx+Sy(s)*ny+Sz(s)*nz)                                     &
    /sqrt(Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))
    shade=abs(shade)
    shade=0.4+0.6*shade

    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if(c2 < 0 .or. &
      ( abs(Dx(s))+abs(Dy(s))+abs(Dz(s)) ) > 0. ) then 

      x1 = grid % nodes(SideN(s,1)) % x
      x2 = grid % nodes(SideN(s,2)) % x
      x3 = grid % nodes(SideN(s,3)) % x
      x4 = grid % nodes(SideN(s,4)) % x

      y1 = grid % nodes(SideN(s,1)) % y
      y2 = grid % nodes(SideN(s,2)) % y
      y3 = grid % nodes(SideN(s,3)) % y
      y4 = grid % nodes(SideN(s,4)) % y

      z1 = grid % nodes(SideN(s,1)) % z
      z2 = grid % nodes(SideN(s,2)) % z
      z3 = grid % nodes(SideN(s,3)) % z
      z4 = grid % nodes(SideN(s,4)) % z

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

      if(SideN(s,0) == 4)                                             &
      write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A7,3F5.2,A15)')       &
                'gs np 0 slw ',                                       &
                int(sclf*xp1),int(sclf*yp1), ' m ',                   & 
                int(sclf*xp2),int(sclf*yp2), ' l ',                   &
                int(sclf*xp3),int(sclf*yp3), ' l ',                   &
                int(sclf*xp4),int(sclf*yp4), ' l ',                   &
                ' cp gs ',                                            &
                red(proces(NewC(c1))),                                &
                green(proces(NewC(c1))),                              &
                blue(proces(NewC(c1))),                               &
                ' srgb f gr s gr'

      if(SideN(s,0) == 3)                                             &
      write(9,'(A12,2I8,A3,2I8,A3,2I8,A3,A7,3F5.2,A15)')              &
                'gs np 0 slw ',                                       &
                int(sclf*xp1),int(sclf*yp1), ' m ',                   & 
                int(sclf*xp2),int(sclf*yp2), ' l ',                   &
                int(sclf*xp3),int(sclf*yp3), ' l ',                   &
                ' cp gs ',                                            &
                red(proces(NewC(c1))),                                &
                green(proces(NewC(c1))),                              &
                blue(proces(NewC(c1))),                               &
                ' srgb f gr s gr'
    end if

  end do 

  ix1 = xmaxb + boxsize / 2
  iy1 = ymaxb
  ix2 = ix1+boxsize
  iy2 = iy1+boxsize
  do n=1,n_sub
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
    write(9,'(A6,2I8,A12,I3, A1, A6)')                         &
              'gs np ',                                        &
              int(1.0/sclp*(ix2+boxsize/2)),int(1.0/sclp*iy1), &
              ' m (process:', n, ')',                          &
              ' sh gr'               
  end do
  write(9,'(A4)') ' gr '

  close(9)

  goto 1

  deallocate(indx)
  deallocate(work)

  end subroutine Save_Eps_Decomposed
