!==============================================================================!
  subroutine Smooth_Grid
!------------------------------------------------------------------------------!
!   Smooths the grid lines by a Laplacian-like algorythm.                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i, j, k, m
  real                 :: x_new_tmp, y_new_tmp, z_new_tmp    
  real                 :: x_max, y_max, z_max, x_min, y_min, z_min 
  real                 :: x1, y1, z1, x8, y8, z8 
  integer              :: reg
  real, allocatable    :: x_node_new(:), y_node_new(:), z_node_new(:) 
  integer, allocatable :: node_to_nodes(:,:)
!==============================================================================!

  ! Allocate memory for additional arrays
  allocate(x_node_new(MAXN)); x_node_new=0
  allocate(y_node_new(MAXN)); y_node_new=0
  allocate(z_node_new(MAXN)); z_node_new=0
  allocate(node_to_nodes(MAXN,0:40)); node_to_nodes=0

  write(*,*) '# Now smoothing the cells. This may take a while !' 

  do n=1,Nn
    node_to_nodes(n,0) = 0
  end do 

  x_max=-HUGE
  y_max=-HUGE
  z_max=-HUGE
  x_min=+HUGE
  y_min=+HUGE
  z_min=+HUGE
  do n=1,Nn
    x_max=max(x_node(n),x_max) 
    y_max=max(y_node(n),y_max) 
    z_max=max(z_node(n),z_max) 
    x_min=min(x_node(n),x_min) 
    y_min=min(y_node(n),y_min) 
    z_min=min(z_node(n),z_min) 
  end do 

  !-----------------------!
  !   Connect the nodes   !
  !-----------------------!
  do c=1,Nc            ! through cells
    do i=1,8           ! through nodes of a cell 
      n = CellN(c,i)   ! first cell
      do j=1,8         ! through nodes of a cell 
        m = CellN(c,j) ! second cell 
        if(n  /=  m) then 
          do k=1,node_to_nodes(n,0)
            if(node_to_nodes(n,k) == m) goto 10            
          end do
          node_to_nodes(n,0)=node_to_nodes(n,0)+1
          node_to_nodes(n,node_to_nodes(n,0)) = m
        end if
10      end do
    end do
  end do 

  !----------------------------!
  !   Browse through regions   !
  !----------------------------!
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
          node_to_nodes(n,0) = 0
        endif
      end do
    end if
  end do

  !---------------------!
  !   Smooth the grid   !
  !---------------------!
  do reg=1,n_smoothing_regions
    write(*,*) '# Now smoothing region ',reg,' with:',              &
                smooth_iters(reg), ' iterations.'

    do j=1,smooth_iters(reg)         

      ! Calculate new coordinates using the 
      ! old values (x_node(),y_node(),z_node())
      do n=1,Nn
        if(node_to_nodes(n,0)   >  0) then
          x_new_tmp=0.0
          y_new_tmp=0.0
          z_new_tmp=0.0
          do i=1,node_to_nodes(n,0)
            x_new_tmp = x_new_tmp + x_node(node_to_nodes(n,i))
            y_new_tmp = y_new_tmp + y_node(node_to_nodes(n,i))
            z_new_tmp = z_new_tmp + z_node(node_to_nodes(n,i))
          end do
          x_new_tmp = x_new_tmp / (1.0*node_to_nodes(n,0)) 
          y_new_tmp = y_new_tmp / (1.0*node_to_nodes(n,0)) 
          z_new_tmp = z_new_tmp / (1.0*node_to_nodes(n,0)) 
          if(x_node(n)  > 0.001*x_min .and. x_node(n)  < 0.999*x_max)           &
          x_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*x_node(n)                &
                  +      smooth_relax(reg)*walln(n) *x_new_tmp
          if(y_node(n)  > 0.001*y_min .and. y_node(n)  < 0.999*y_max)           &
          y_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*y_node(n)                &
                  +      smooth_relax(reg)*walln(n) *y_new_tmp
          if(z_node(n)  > 0.001*z_min .and. z_node(n)  < 0.999*z_max)           &
          z_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*z_node(n)                &
                  +      smooth_relax(reg)*walln(n)* z_new_tmp
        end if
      end do

      ! Update coordinates
      do n=1,Nn
        if(node_to_nodes(n,0)   >  0) then

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
              if(x_node(n)  > 0.001*x_min .and. x_node(n)  < 0.999*x_max)       &
              x_node(n)=x_node_new(n)
            end if

            if(smooth_in_y(reg)) then
              if(y_node(n)  > 0.001*y_min .and. y_node(n)  < 0.999*y_max)       &
              y_node(n)=y_node_new(n)
            end if 

            if(smooth_in_z(reg)) then
              if(z_node(n)  > 0.001*z_min .and. z_node(n)  < 0.999*z_max)       &
              z_node(n)=z_node_new(n)
            end if 

          end if  ! if the point belongs to region
        end if    ! if the point is not excluded from smoothing
      end do      ! through nodes
    end do
  end do

  deallocate(x_node_new)
  deallocate(y_node_new)
  deallocate(z_node_new)
  deallocate(node_to_nodes)

  end subroutine Smooth_Grid
