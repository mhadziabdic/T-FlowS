!==============================================================================!
  subroutine Save_Eps_Cut(grid, face_g_dx, face_g_dy, dir)
!------------------------------------------------------------------------------!
!   Writes grid in encapsulated postscript format.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use all_mod, only: name
  use gen_mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: xg(grid % max_n_nodes),  &  ! if I recall correctly (BN)
                     yg(grid % max_n_nodes),  &  ! these were supposed to be
                     zg(grid % max_n_nodes)      ! general coordnates
  real            :: face_g_dx(grid % max_n_faces),  &
                     face_g_dy(grid % max_n_faces) 
  character       :: dir
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c1, c2, count, lw
  character(len=80) :: name_eps, answer
  real              :: sclf, sclp, xmax,xmin,ymax,ymin,zmax,zmin, z0, fc
  real              :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xin(4),yin(4) 
!==============================================================================!

  write(*,'(A39)')        ' #==============================================='
  write(*,'(A10,A1,A28)') ' # Making ', dir, '.eps cut through the domain'
  write(*,'(A39)')        ' #-----------------------------------------------'

  write(*,*) '# Enter the ',dir,' coordinate for cutting or type skip to exit: '
  call Tokenizer_Mod_Read_Line(5)
  read(line % whole, *) answer
  call To_Upper_Case(answer)
  if(answer == 'SKIP') return  
  read(answer, *) z0
  if(z0 == 0) return
  write(*,*) '# Z0 = ', z0

  !----------------------!
  !                      !
  !   Create .eps file   !
  !                      !
  !----------------------!
  name_eps = name
  name_eps(len_trim(name)+1:len_trim(name)+6) = '. .eps'
  name_eps(len_trim(name)+2:len_trim(name)+2) = dir
  write(*, *) '# Now creating the file:', name_eps

  if(dir == 'x') then
    xmax=maxval(grid % yn(1:grid % n_nodes))
    ymax=maxval(grid % zn(1:grid % n_nodes))
    zmax=maxval(grid % xn(1:grid % n_nodes))
    xmin=minval(grid % yn(1:grid % n_nodes))
    ymin=minval(grid % zn(1:grid % n_nodes))
    zmin=minval(grid % xn(1:grid % n_nodes))
  else if(dir == 'y') then
    xmax=maxval(grid % zn(1:grid % n_nodes))
    ymax=maxval(grid % xn(1:grid % n_nodes))
    zmax=maxval(grid % yn(1:grid % n_nodes))
    xmin=minval(grid % zn(1:grid % n_nodes))
    ymin=minval(grid % xn(1:grid % n_nodes))
    zmin=minval(grid % yn(1:grid % n_nodes))
  else if(dir == 'z') then
    xmax=maxval(grid % xn(1:grid % n_nodes))
    ymax=maxval(grid % yn(1:grid % n_nodes))
    zmax=maxval(grid % zn(1:grid % n_nodes))
    xmin=minval(grid % xn(1:grid % n_nodes))
    ymin=minval(grid % yn(1:grid % n_nodes))
    zmin=minval(grid % zn(1:grid % n_nodes))
  end if

  sclf = 100000.0/max((xmax-xmin),(ymax-ymin))
  sclp = 0.005 

  open(9, file=name_eps)

  write(9, '(A24)') '%!PS-Adobe-2.0 EPSF-1.2 '
  write(9, '(A24)') '%% Created by:  TFlowS %%'
  write(9, '(A15,4I7)') '%%BoundingBox: ',                          &
    int(xmin*sclf*sclp)-1,int(ymin*sclf*sclp)-1,                    &
    int(xmax*sclf*sclp)+1,int(ymax*sclf*sclp)+1 
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
  write(9, '(A24)') '/slw  {setlinewidth} def'
  write(9, '(A24)') '% End of brevations.    '

  write(9,'(A4)') ' gs '
  write(9, '(A24)') '% Scale:                '
  write(9, '(2F10.6,A4)') sclp,sclp, ' sc '

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    x1 = xg(grid % faces_n(1,s))
    y1 = yg(grid % faces_n(1,s))
    z1 = zg(grid % faces_n(1,s))

    x2 = xg(grid % faces_n(2,s))
    y2 = yg(grid % faces_n(2,s))
    z2 = zg(grid % faces_n(2,s))

    x3 = xg(grid % faces_n(3,s))
    y3 = yg(grid % faces_n(3,s))
    z3 = zg(grid % faces_n(3,s))

    x4 = xg(grid % faces_n(4,s))
    y4 = yg(grid % faces_n(4,s))
    z4 = zg(grid % faces_n(4,s))

    count=0
    ! 1-2
    if(z1<=z0.and.z2>=z0 .or. z1>=z0 .and. z2<=z0) then 
      count=count+1
      fc = (z2-z0)/(z2-z1+TINY)
      xin(count) = fc*x1+(1.0-fc)*x2
      yin(count) = fc*y1+(1.0-fc)*y2
    end if
    ! 1-3
    if(z1<=z0 .and. z3>=z0 .or. z1>=z0 .and. z3<=z0) then 
      count=count+1
      fc = (z3-z0)/(z3-z1+TINY)
      xin(count) = fc*x1+(1.0-fc)*x3
      yin(count) = fc*y1+(1.0-fc)*y3
    end if
    ! 2-3
    if(z2<=z0 .and. z3>=z0 .or. z2>=z0 .and. z3<=z0) then 
      count=count+1
      fc = (z3-z0)/(z3-z2+TINY)
      xin(count) = fc*x2+(1.0-fc)*x3
      yin(count) = fc*y2+(1.0-fc)*y3
    end if
    ! For quadrilateral faces
    if(grid % faces_n_nodes(s)==4) then
      ! 1-4
      if(z1<=z0 .and. z4>=z0 .or. z1>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z1+TINY)
        xin(count) = fc*x1+(1.0-fc)*x4
        yin(count) = fc*y1+(1.0-fc)*y4
      end if
      ! 2-4
      if(z2<=z0 .and. z4>=z0 .or. z2>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z2+TINY)
        xin(count) = fc*x2+(1.0-fc)*x4
        yin(count) = fc*y2+(1.0-fc)*y4
      end if
      ! 3-4
      if(z3<=z0 .and. z4>=z0 .or. z3>=z0 .and. z4<=z0) then 
        count=count+1
        fc = (z4-z0)/(z4-z3+TINY)
        xin(count) = fc*x3+(1.0-fc)*x4
        yin(count) = fc*y3+(1.0-fc)*y4
      end if
    end if

    ! Line width 
    lw = 0 

    ! If you want to move the domain, do something like:
    ! xin = (1.0 - xin) 

    if(count == 2) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,A8)')                        &
                'gs np ',lw,' slw ',                                &
                int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
                int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
                ' cp s gr' 
      if( grid % dx(s) /= 0.0 .or. grid % dy(s) /= 0.0 ) then
        write(9,'(A6,I2,A5,2I8,A3,2I8,A3,A8)')                      &
                  'gs np ',lw,' slw ',                              &
                  int(sclf*(xin(1)-face_g_dx(s))),                  &
                  int(sclf*(yin(1)-face_g_dy(s))), ' m ',           &
                  int(sclf*(xin(2)-face_g_dx(s))),                  &
                  int(sclf*(yin(2)-face_g_dy(s))), ' l ',           &
                  ' cp s gr' 
      end if
    end if

    if(count == 3) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,2I8,A3,A8)')                 &
                'gs np ',lw,' slw ',                                &
                int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
                int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
                int(sclf*xin(3)),int(sclf*yin(3)), ' l ',           &
                ' cp s gr' 
    end if

    if(count == 4) then
      write(9,'(A6,I2,A5,2I8,A3,2I8,A3,2I8,A3,2I8,A3,A8)')          &
                'gs np ',lw,' slw ',                                &
                int(sclf*xin(1)),int(sclf*yin(1)), ' m ',           &
                int(sclf*xin(2)),int(sclf*yin(2)), ' l ',           &
                int(sclf*xin(3)),int(sclf*yin(3)), ' l ',           &
                int(sclf*xin(4)),int(sclf*yin(4)), ' l ',           &
                ' cp s gr' 
    end if

  end do 

  write(9,'(A4)') ' gr '

  close(9)

  end subroutine
