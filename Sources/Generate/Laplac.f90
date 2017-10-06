!======================================================================!
  subroutine Laplac(b,i,j,k,wx16,wx24,wx35,                         &
                            wy16,wy24,wy35,                         &
                            wz16,wz24,wz35)
!----------------------------------------------------------------------!
!   Places the nodes inside the block using Laplace-like function      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------! 
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: b, i, j, k
  real    :: wx16, wx24, wx35, wy16, wy24, wy35, wz16, wz24, wz35 
!-------------------------------[Locals]-------------------------------!
  real    :: xt(8), yt(8), zt(8)
  integer :: NI, NJ, NK, n, n1, n2, n3, n4, n5, n6
  real    :: xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3 
  real    :: xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6
!======================================================================!

  NI=block_resolutions(b,1)
  NJ=block_resolutions(b,2)
  NK=block_resolutions(b,3)   

  do n=1,8
    xt(n) = x_point(block_points(b, n)) 
    yt(n) = y_point(block_points(b, n))
    zt(n) = z_point(block_points(b, n))
  end do  

  n = NN + (k-1)*NI*NJ + (j-1)*NI + i

!----- node numbers at the block faces
  n1 = NN + ( 1-1)*NI*NJ + ( j-1)*NI + i     !  ->  k == 1
  n2 = NN + ( k-1)*NI*NJ + ( j-1)*NI + 1     !  ->  i == 1
  n3 = NN + ( k-1)*NI*NJ + ( 1-1)*NI + i     !  ->  j == 1
  n4 = NN + ( k-1)*NI*NJ + ( j-1)*NI + NI    !  ->  i == NI
  n5 = NN + ( k-1)*NI*NJ + (NJ-1)*NI + i     !  ->  j == NJ
  n6 = NN + (NK-1)*NI*NJ + ( j-1)*NI + i     !  ->  k == NK

!----- face I
  if(x_node(n1) == HUGE) then
    xf1=( ((NI-i)*xt(1) + (i-1)*xt(2)) * (NJ-j) +                   &
          ((NI-i)*xt(3) + (i-1)*xt(4)) * (j-1)  ) /((NI-1)*(NJ-1))
    yf1=( ((NI-i)*yt(1) + (i-1)*yt(2)) * (NJ-j) +                   &
          ((NI-i)*yt(3) + (i-1)*yt(4)) * (j-1)  ) /((NI-1)*(NJ-1))
    zf1=( ((NI-i)*zt(1) + (i-1)*zt(2)) * (NJ-j) +                   &
          ((NI-i)*zt(3) + (i-1)*zt(4)) * (j-1)  ) /((NI-1)*(NJ-1))
  else
    xf1=x_node(n1)
    yf1=y_node(n1)
    zf1=z_node(n1)
  end if

!----- face VI
  if(x_node(n6) == HUGE) then
    xf6=( ((NI-i)*xt(5) + (i-1)*xt(6)) * (NJ-j) +                   &
          ((NI-i)*xt(7) + (i-1)*xt(8)) * (j-1)  ) /((NI-1)*(NJ-1))
    yf6=( ((NI-i)*yt(5) + (i-1)*yt(6)) * (NJ-j) +                   &
          ((NI-i)*yt(7) + (i-1)*yt(8)) * (j-1)  ) /((NI-1)*(NJ-1))
    zf6=( ((NI-i)*zt(5) + (i-1)*zt(6)) * (NJ-j) +                   &
          ((NI-i)*zt(7) + (i-1)*zt(8)) * (j-1)  ) /((NI-1)*(NJ-1))
  else
    xf6=x_node(n6)
    yf6=y_node(n6)
    zf6=z_node(n6)
  endif

!----- face III
  if(x_node(n3) == HUGE) then
    xf3=( ((NI-i)*xt(1) + (i-1)*xt(2)) * (Nk-k) +                   &
          ((NI-i)*xt(5) + (i-1)*xt(6)) * (k-1)  ) /((NI-1)*(NK-1))
    yf3=( ((NI-i)*yt(1) + (i-1)*yt(2)) * (Nk-k) +                   &
          ((NI-i)*yt(5) + (i-1)*yt(6)) * (k-1)  ) /((NI-1)*(NK-1))
    zf3=( ((NI-i)*zt(1) + (i-1)*zt(2)) * (Nk-k) +                   &
          ((NI-i)*zt(5) + (i-1)*zt(6)) * (k-1)  ) /((NI-1)*(NK-1))
  else
    xf3=x_node(n3)
    yf3=y_node(n3)
    zf3=z_node(n3)
  endif

!----- face V
  if(x_node(n5) == HUGE) then
    xf5=( ((NI-i)*xt(3) + (i-1)*xt(4)) * (Nk-k) +                   &
          ((NI-i)*xt(7) + (i-1)*xt(8)) * (k-1)  ) /((NI-1)*(NK-1))
    yf5=( ((NI-i)*yt(3) + (i-1)*yt(4)) * (Nk-k) +                   &
          ((NI-i)*yt(7) + (i-1)*yt(8)) * (k-1)  ) /((NI-1)*(NK-1))
    zf5=( ((NI-i)*zt(3) + (i-1)*zt(4)) * (Nk-k) +                   &
          ((NI-i)*zt(7) + (i-1)*zt(8)) * (k-1)  ) /((NI-1)*(NK-1))
  else
    xf5=x_node(n5)
    yf5=y_node(n5)
    zf5=z_node(n5)
  endif

!----- face II
  if(x_node(n2) == HUGE) then
    xf2=( ((NJ-j)*xt(1) + (j-1)*xt(3)) * (NK-k) +                   &
          ((NJ-j)*xt(5) + (j-1)*xt(7)) * (k-1)  ) /((NJ-1)*(NK-1))
    yf2=( ((NJ-j)*yt(1) + (j-1)*yt(3)) * (NK-k) +                   &
          ((NJ-j)*yt(5) + (j-1)*yt(7)) * (k-1)  ) /((NJ-1)*(NK-1))
    zf2=( ((NJ-j)*zt(1) + (j-1)*zt(3)) * (NK-k) +                   &
          ((NJ-j)*zt(5) + (j-1)*zt(7)) * (k-1)  ) /((NJ-1)*(NK-1))
  else
    xf2=x_node(n2)
    yf2=y_node(n2)
    zf2=z_node(n2)
  endif

!----- face IV
  if(x_node(n4) == HUGE) then
    xf4=( ((NJ-j)*xt(2) + (j-1)*xt(4)) * (NK-k) +                   &
          ((NJ-j)*xt(6) + (j-1)*xt(8)) * (k-1)  ) /((NJ-1)*(NK-1))
    yf4=( ((NJ-j)*yt(2) + (j-1)*yt(4)) * (NK-k) +                   &
          ((NJ-j)*yt(6) + (j-1)*yt(8)) * (k-1)  ) /((NJ-1)*(NK-1))
    zf4=( ((NJ-j)*zt(2) + (j-1)*zt(4)) * (NK-k) +                   &
          ((NJ-j)*zt(6) + (j-1)*zt(8)) * (k-1)  ) /((NJ-1)*(NK-1))
  else
    xf4=x_node(n4)
    yf4=y_node(n4)
    zf4=z_node(n4)
  endif

  if( x_node(n) == HUGE ) then
    x_node(n) = ( xf1*(NK-k) + xf6*(k-1) ) * wx16 / (NK-1)               &
              + ( xf2*(NI-i) + xf4*(i-1) ) * wx24 / (NI-1)               &
              + ( xf3*(NJ-j) + xf5*(j-1) ) * wx35 / (NJ-1) 

    y_node(n) = ( yf1*(NK-k) + yf6*(k-1) ) * wy16 / (NK-1)               &
              + ( yf2*(NI-i) + yf4*(i-1) ) * wy24 / (NI-1)               &
              + ( yf3*(NJ-j) + yf5*(j-1) ) * wy35 / (NJ-1) 

    z_node(n) = ( zf1*(NK-k) + zf6*(k-1) ) * wz16 / (NK-1)               &
              + ( zf2*(NI-i) + zf4*(i-1) ) * wz24 / (NI-1)               &
              + ( zf3*(NJ-j) + zf5*(j-1) ) * wz35 / (NJ-1) 
  end if

  end subroutine Laplac
