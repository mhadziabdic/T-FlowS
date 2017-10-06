!======================================================================!
  subroutine Linija(b, w, is,js,ks,ie,je,ke)
!----------------------------------------------------------------------!
!   Places the nodes on the line defined with local block position     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: b, is,js,ks,ie,je,ke
  real    :: w
!------------------------------[Calling]-------------------------------!
  real :: atanh
!-------------------------------[Locals]-------------------------------!
  integer :: N, ni, nj, nk, i, j, k, node, case
  real    :: x0, y0, z0, delx, dely, delz, t, dt, ddt, pr, xi
!======================================================================!

  ni=block_resolutions(b,1)
  nj=block_resolutions(b,2)
  nk=block_resolutions(b,3)   

  x0   = x_node(NN+(ks-1)*ni*nj+(js-1)*ni+is)
  y0   = y_node(NN+(ks-1)*ni*nj+(js-1)*ni+is)
  z0   = z_node(NN+(ks-1)*ni*nj+(js-1)*ni+is)
  delx = x_node(NN+(ke-1)*ni*nj+(je-1)*ni+ie)-x0
  dely = y_node(NN+(ke-1)*ni*nj+(je-1)*ni+ie)-y0
  delz = z_node(NN+(ke-1)*ni*nj+(je-1)*ni+ie)-z0

  N=max( (ie-is), (je-js),  (ke-ks) )

!-----------------------------!
!     Linear distribution     !
!-----------------------------!
  if(w  > 0.0) then
    ddt = ( 2.0*(1.0-w) ) / ( 1.0*N*(1.0*N-1.0)*(1.0+w) )
    t=0.0
    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie /= is ) then
            dt=1.0/(1.0*N)+(1.0*i-0.5*(1.0*N+1)) * ddt
            t=t+dt
            node = NN + (k-1)*ni*nj + (j-1)*ni + i+1
            if( (i  < ie).and.(x_node(node) == HUGE) ) then 
              x_node(node) = x0 + t*delx
              y_node(node) = y0 + t*dely
              z_node(node) = z0 + t*delz
            endif
          end if 
          if( je /= js ) then
            dt=1.0/(1.0*N)+(1.0*j-0.5*(1.0*N+1)) * ddt
            t=t+dt
            node = NN + (k-1)*ni*nj + (j-0)*ni + i 
            if( (j  < je).and.(x_node(node) == HUGE) ) then 
              x_node(node) = x0 + t*delx
              y_node(node) = y0 + t*dely
              z_node(node) = z0 + t*delz
            endif
          end if 
          if( ke /= ks ) then
            dt=1.0/(1.0*N)+(1.0*k-0.5*(1.0*N+1)) * ddt
            t=t+dt
            node = NN + (k-0)*ni*nj + (j-1)*ni + i 
            if( (k  < ke).and.(x_node(node) == HUGE) ) then 
              x_node(node) = x0 + t*delx
              y_node(node) = y0 + t*dely
              z_node(node) = z0 + t*delz
            endif
          end if 
        end do
      end do
    end do

!---------------------------------!
!     Hyperbolic distribution     !
!---------------------------------!
  else
    case = 0
    if     ((w  >  -0.5).and.(w <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(w))    
      case = 1
    else if((w >=  -0.75).and.(w  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(w))            
      case = 2
    else
      pr = -w
      case = 3 
    endif

    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie /= is ) then
    if(case == 1) xi = -1.0*(1.0*i)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*i)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*i)/(1.0*N)
            node = NN + (k-1)*ni*nj + (j-1)*ni + i+1
            if( (i  < ie).and.(x_node(node) == HUGE) ) then 
    if    (case == 1) then
      x_node(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x_node(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x_node(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
            endif
          end if 
          if( je /= js ) then
    if(case == 1) xi = -1.0*(1.0*j)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*j)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*j)/(1.0*N)
            node = NN + (k-1)*ni*nj + (j-0)*ni + i 
            if( (j  < je).and.(x_node(node) == HUGE) ) then 
    if    (case == 1) then
      x_node(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x_node(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x_node(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
            endif
          end if 
          if( ke /= ks ) then
    if(case == 1) xi = -1.0*(1.0*k)/(1.0*N)
    if(case == 2) xi =  1.0 - 1.0*(1.0*k)/(1.0*N)
    if(case == 3) xi = -1.0 + 2.0*(1.0*k)/(1.0*N)
            node = NN + (k-0)*ni*nj + (j-1)*ni + i 
            if( (k  < ke).and.(x_node(node) == HUGE) ) then 
    if    (case == 1) then
      x_node(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 2) then
      x_node(node) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + dely - (tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + delz - (tanh(xi*atanh(pr))/pr)*delz
    elseif(case == 3) then
      x_node(node) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      y_node(node) = y0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
      z_node(node) = z0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
    endif 
            endif
          end if 
        end do
      end do
    end do

  endif

  end subroutine Linija
