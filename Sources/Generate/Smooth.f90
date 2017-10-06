!======================================================================!
  subroutine Smooth
!----------------------------------------------------------------------!
!   Smooths the grid lines by a Laplacian-like algorythm.              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer             :: c, n, i,j,k,m
  real                :: xTnew, yTnew, zTnew    ! temporary new values
  real                :: xmax, ymax, zmax, xmin, ymin, zmin 
  real                :: x1,y1,z1,x8,y8,z8 
  integer             :: reg
  real,allocatable    :: xnew(:), ynew(:), znew(:) ! new values
  integer,allocatable :: NodeN(:,:)
!======================================================================!

!---- allocate memory for additional arrays
  allocate(xnew(MAXN)); xnew=0
  allocate(ynew(MAXN)); ynew=0
  allocate(znew(MAXN)); znew=0
  allocate(NodeN(MAXN,0:40)); NodeN=0

  write(*,*) 'Now smoothing the cells. This may take a while !' 

  do n=1,Nn
    NodeN(n,0) = 0
  end do 

  xmax=-HUGE
  ymax=-HUGE
  zmax=-HUGE
  xmin=+HUGE
  ymin=+HUGE
  zmin=+HUGE
  do n=1,Nn
    xmax=max(x_node(n),xmax) 
    ymax=max(y_node(n),ymax) 
    zmax=max(z_node(n),zmax) 
    xmin=min(x_node(n),xmin) 
    ymin=min(y_node(n),ymin) 
    zmin=min(z_node(n),zmin) 
  end do 

!---------------------------!
!     Connect the nodes     !
!---------------------------!
  do c=1,Nc            ! through cells
    do i=1,8           ! through nodes of a cell 
      n = CellN(c,i)   ! first cell
      do j=1,8         ! through nodes of a cell 
        m = CellN(c,j) ! second cell 
        if(n  /=  m) then 
          do k=1,NodeN(n,0)
            if(NodeN(n,k) == m) goto 10            
          end do
          NodeN(n,0)=NodeN(n,0)+1
          NodeN(n,NodeN(n,0)) = m
        end if
10      end do
    end do
  end do 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Browse through regions     !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  do reg=1,n_smoothing_regions
    if( ( .not. smooth_in_x(reg) ) .and.  &
        ( .not. smooth_in_y(reg) ) .and.  &
        ( .not. smooth_in_z(reg) ) ) then
      do n=1,Nn
        x1=smooth_regions(reg,1)
        y1=smooth_regions(reg,2)
        z1=smooth_regions(reg,3)
        x8=smooth_regions(reg,4)
        y8=smooth_regions(reg,5)
        z8=smooth_regions(reg,6)
        if( (x1 <= x_node(n)) .and. (x_node(n) <= x8) .and.  &
            (y1 <= y_node(n)) .and. (y_node(n) <= y8) .and.  &
            (z1 <= z_node(n)) .and. (z_node(n) <= z8) ) then
          NodeN(n,0) = 0
        endif
      end do
    end if
  end do

!->>> do n=1,Nn
!->>>   write(*,*) '=>', n, NodeN(n,0)
!->>>   write(*,*) (NodeN(n,i), i=1,NodeN(n,0) )
!->>> end do

!-------------------------!
!     Smooth the grid     !
!-------------------------!
  do reg=1,n_smoothing_regions
    write(*,*) 'Now smoothing region ',reg,' with:',              &
                smooth_iters(reg), ' iterations.'

    do j=1,smooth_iters(reg)         

!---- calculate new coordinates using the old values (x_node(),y_node(),z_node())
      do n=1,Nn
        if(NodeN(n,0)   >  0) then
          xTnew=0.0
          yTnew=0.0
          zTnew=0.0
          do i=1,NodeN(n,0)
            xTnew = xTnew + x_node(NodeN(n,i))
            yTnew = yTnew + y_node(NodeN(n,i))
            zTnew = zTnew + z_node(NodeN(n,i))
          end do
          xTnew = xTnew / (1.0*NodeN(n,0)) 
          yTnew = yTnew / (1.0*NodeN(n,0)) 
          zTnew = zTnew / (1.0*NodeN(n,0)) 
          if(x_node(n)  > 0.001*xmin .and. x_node(n)  < 0.999*xmax)           &
          xnew(n) = (1.0-smooth_relax(reg)*walln(n))*x_node(n)                &
                  +      smooth_relax(reg)*walln(n) *xTnew
          if(y_node(n)  > 0.001*ymin .and. y_node(n)  < 0.999*ymax)           &
          ynew(n) = (1.0-smooth_relax(reg)*walln(n))*y_node(n)                &
                  +      smooth_relax(reg)*walln(n) *yTnew
          if(z_node(n)  > 0.001*zmin .and. z_node(n)  < 0.999*zmax)           &
          znew(n) = (1.0-smooth_relax(reg)*walln(n))*z_node(n)                &
                  +      smooth_relax(reg)*walln(n)* zTnew
        end if
      end do

!---- update coordinates
      do n=1,Nn
        if(NodeN(n,0)   >  0) then

          x1=smooth_regions(reg,1)
          y1=smooth_regions(reg,2)
          z1=smooth_regions(reg,3)
          x8=smooth_regions(reg,4)
          y8=smooth_regions(reg,5)
          z8=smooth_regions(reg,6)

          if( (x1 <= x_node(n)) .and. (x_node(n) <= x8) .and.                 &
              (y1 <= y_node(n)) .and. (y_node(n) <= y8) .and.                 &
              (z1 <= z_node(n)) .and. (z_node(n) <= z8) ) then

            if(smooth_in_x(reg)) then
              if(x_node(n)  > 0.001*xmin .and. x_node(n)  < 0.999*xmax)       &
              x_node(n)=xnew(n)
            end if

            if(smooth_in_y(reg)) then
              if(y_node(n)  > 0.001*ymin .and. y_node(n)  < 0.999*ymax)       &
              y_node(n)=ynew(n)
            end if 

            if(smooth_in_z(reg)) then
              if(z_node(n)  > 0.001*zmin .and. z_node(n)  < 0.999*zmax)       &
              z_node(n)=znew(n)
            end if 

          end if  ! if the point belongs to region
        end if    ! if the point is not excluded from smoothing
      end do      ! through nodes
    end do
  end do

  deallocate(xnew)
  deallocate(ynew)
  deallocate(znew)
  deallocate(NodeN)

  end subroutine Smooth
