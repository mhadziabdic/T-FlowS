!==============================================================================!
  subroutine Save_Eps_Decomposed(grid)
!------------------------------------------------------------------------------!
!   Saves the whole grid in encapsulated postscript.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use div_mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n, s, s0, c1, c2
  integer             :: xmaxb, xminb, ymaxb, yminb, xlegend
  character(len=80)   :: name_eps, answer
  real                :: sclf, sclp, xmax,xmin,ymax,ymin,zmax,zmin
  real                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,        &
                         xk,yk,zk,alfa,beta,gama,nx,ny,nz,shade,fs,  &
                         xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4 
  real                :: red(18), green(18), blue(18)
  integer,allocatable :: indx(:)
  real,allocatable    :: work(:) 
  integer             :: ix1, ix2, iy1, iy2, boxsize
!==============================================================================!

  if(n_sub > 18) then
    print *, '# Error message from Save_Eps_Decomposed: the subroutine '
    print *, '# can''t plot more than 18 partitions.  Skipping it!'
  end if

  ! Allocate the memory
  allocate(indx(grid % n_faces)); indx=0
  allocate(work(grid % n_faces)); work=0   

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
1 print *, '# Enter the camera coordinates (skip to exit): '
  call Tokenizer_Mod_Read_Line(5)
  if(line % n_tokens == 1) then
    read(line % whole, *) answer
    call To_Upper_Case(answer)
    if(answer == 'SKIP') return 
  else if(line % n_tokens == 3) then
    read(line % whole, *) xk, yk, zk 
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
  print *, '# Enter the file name (without extension): '
  call Tokenizer_Mod_Read_Line(5)
  read(line % whole, *) name_eps 
  name_eps(len_trim(name_eps)+1:len_trim(name_eps)+4) = '.eps'
  print *, '# Now creating the file:', trim(name_eps)

  xmax=maxval(grid % xn(1:grid % n_nodes))
  ymax=maxval(grid % yn(1:grid % n_nodes))
  zmax=maxval(grid % zn(1:grid % n_nodes))
  xmin=minval(grid % xn(1:grid % n_nodes))
  ymin=minval(grid % yn(1:grid % n_nodes))
  zmin=minval(grid % zn(1:grid % n_nodes))
  sclf = 100000.0/max((xmax-xmin),(ymax-ymin),(zmax-zmin))
  sclp = 0.005 

  ! Find the bounding box
  xminb=+1000000
  xmaxb=-1000000
  yminb=+1000000
  ymaxb=-1000000
  do n = 1, grid % n_nodes
    if(xk  < 0.0 .and. yk  > 0.0) then
      xp1= - grid % xn(n) * sin(alfa) - grid % yn(n) * sin(beta)
    else if(xk  > 0.0 .and. yk  < 0.0) then
      xp1=   grid % xn(n) * sin(alfa) + grid % yn(n) * sin(beta)
    else if(xk  > 0.0 .and. yk  > 0.0) then
      xp1= - grid % xn(n) * sin(alfa) + grid % yn(n) * sin(beta)
    else
      xp1=   grid % xn(n) * sin(alfa) - grid % yn(n) * sin(beta)
    end if
    xp1=xp1*sclf*sclp
    xmaxb=max(xmaxb,int(xp1))
    xminb=min(xminb,int(xp1))

    yp1=( - grid % xn(n) * cos(alfa)                  &
          - grid % yn(n) * cos(beta) ) * cos(gama) +  &
            grid % zn(n) * sin(gama) 
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

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    indx(s) = s
    work(s) = Distance(xk,yk,zk,                                            &
                fs*grid % xc(c1)+(1.-fs)*(grid % xc(c2)+grid % dx(s)),  &
                fs*grid % yc(c1)+(1.-fs)*(grid % yc(c2)+grid % dy(s)),  &
                fs*grid % zc(c1)+(1.-fs)*(grid % zc(c2)+grid % dz(s)) )
  end do
  call Sort_Real_Carry_Int(work,indx,grid % n_faces,-2)

  do s0 = 1, grid % n_faces
    s=indx(s0)

    shade = (grid % sx(s)*nx + grid % sy(s)*ny + grid % sz(s)*nz)  &
          / sqrt(  grid % sx(s)*grid % sx(s)                       &
                 + grid % sy(s)*grid % sy(s)                       &
                 + grid % sz(s)*grid % sz(s)  )
    shade=abs(shade)
    shade=0.4+0.6*shade

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0 .or. &
      ( abs(grid % dx(s))+abs(grid % dy(s))+abs(grid % dz(s)) ) > 0. ) then 

      x1 = grid % xn(grid % faces_n(1,s))
      x2 = grid % xn(grid % faces_n(2,s))
      x3 = grid % xn(grid % faces_n(3,s))
      x4 = grid % xn(grid % faces_n(4,s))

      y1 = grid % yn(grid % faces_n(1,s))
      y2 = grid % yn(grid % faces_n(2,s))
      y3 = grid % yn(grid % faces_n(3,s))
      y4 = grid % yn(grid % faces_n(4,s))

      z1 = grid % zn(grid % faces_n(1,s))
      z2 = grid % zn(grid % faces_n(2,s))
      z3 = grid % zn(grid % faces_n(3,s))
      z4 = grid % zn(grid % faces_n(4,s))

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

      if(grid % faces_n_nodes(s) == 4)                                &
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

      if(grid % faces_n_nodes(s) == 3)                                &
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

  end subroutine
