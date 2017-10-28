!==============================================================================!
  subroutine Laplac(b,i,j,k,wx16,wx24,wx35,                         &
                            wy16,wy24,wy35,                         &
                            wz16,wz24,wz35)
!------------------------------------------------------------------------------!
!   Places the nodes inside the block using Laplace-like function              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: b, i, j, k
  real    :: wx16, wx24, wx35, wy16, wy24, wy35, wz16, wz24, wz35 
!-----------------------------------[Locals]-----------------------------------!
  real    :: xt(8), yt(8), zt(8)
  integer :: ni, nj, nk, n, n1, n2, n3, n4, n5, n6
  real    :: xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3 
  real    :: xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6
!==============================================================================!

  ni = dom % blocks(b) % resolutions(1)
  nj = dom % blocks(b) % resolutions(2)
  nk = dom % blocks(b) % resolutions(3)   

  do n=1,8
    xt(n) = dom % points( dom % blocks(b) % corners(n) ) % x
    yt(n) = dom % points( dom % blocks(b) % corners(n) ) % y
    zt(n) = dom % points( dom % blocks(b) % corners(n) ) % z
  end do  

  n = NN + (k-1)*ni*nj + (j-1)*ni + i

  ! Node numbers at the block faces
  n1 = NN + ( 1-1)*ni*nj + ( j-1)*ni + i     !  ->  k == 1
  n2 = NN + ( k-1)*ni*nj + ( j-1)*ni + 1     !  ->  i == 1
  n3 = NN + ( k-1)*ni*nj + ( 1-1)*ni + i     !  ->  j == 1
  n4 = NN + ( k-1)*ni*nj + ( j-1)*ni + ni    !  ->  i == ni
  n5 = NN + ( k-1)*ni*nj + (nj-1)*ni + i     !  ->  j == nj
  n6 = NN + (nk-1)*ni*nj + ( j-1)*ni + i     !  ->  k == nk

  ! Face I
  if(grid % nodes(n1) % x == HUGE) then
    xf1=( ((ni-i)*xt(1) + (i-1)*xt(2)) * (nj-j) +                   &
          ((ni-i)*xt(3) + (i-1)*xt(4)) * (j-1)  ) /((ni-1)*(nj-1))
    yf1=( ((ni-i)*yt(1) + (i-1)*yt(2)) * (nj-j) +                   &
          ((ni-i)*yt(3) + (i-1)*yt(4)) * (j-1)  ) /((ni-1)*(nj-1))
    zf1=( ((ni-i)*zt(1) + (i-1)*zt(2)) * (nj-j) +                   &
          ((ni-i)*zt(3) + (i-1)*zt(4)) * (j-1)  ) /((ni-1)*(nj-1))
  else
    xf1 = grid % nodes(n1) % x
    yf1 = grid % nodes(n1) % y
    zf1 = grid % nodes(n1) % z
  end if

  ! Face VI
  if(grid % nodes(n6) % x == HUGE) then
    xf6=( ((ni-i)*xt(5) + (i-1)*xt(6)) * (nj-j) +                   &
          ((ni-i)*xt(7) + (i-1)*xt(8)) * (j-1)  ) /((ni-1)*(nj-1))
    yf6=( ((ni-i)*yt(5) + (i-1)*yt(6)) * (nj-j) +                   &
          ((ni-i)*yt(7) + (i-1)*yt(8)) * (j-1)  ) /((ni-1)*(nj-1))
    zf6=( ((ni-i)*zt(5) + (i-1)*zt(6)) * (nj-j) +                   &
          ((ni-i)*zt(7) + (i-1)*zt(8)) * (j-1)  ) /((ni-1)*(nj-1))
  else
    xf6 = grid % nodes(n6) % x
    yf6 = grid % nodes(n6) % y
    zf6 = grid % nodes(n6) % z
  endif

  ! Face III
  if(grid % nodes(n3) % x == HUGE) then
    xf3=( ((ni-i)*xt(1) + (i-1)*xt(2)) * (Nk-k) +                   &
          ((ni-i)*xt(5) + (i-1)*xt(6)) * (k-1)  ) /((ni-1)*(nk-1))
    yf3=( ((ni-i)*yt(1) + (i-1)*yt(2)) * (Nk-k) +                   &
          ((ni-i)*yt(5) + (i-1)*yt(6)) * (k-1)  ) /((ni-1)*(nk-1))
    zf3=( ((ni-i)*zt(1) + (i-1)*zt(2)) * (Nk-k) +                   &
          ((ni-i)*zt(5) + (i-1)*zt(6)) * (k-1)  ) /((ni-1)*(nk-1))
  else
    xf3 = grid % nodes(n3) % x
    yf3 = grid % nodes(n3) % y
    zf3 = grid % nodes(n3) % z
  endif

  ! Face V
  if(grid % nodes(n5) % x == HUGE) then
    xf5=( ((ni-i)*xt(3) + (i-1)*xt(4)) * (Nk-k) +                   &
          ((ni-i)*xt(7) + (i-1)*xt(8)) * (k-1)  ) /((ni-1)*(nk-1))
    yf5=( ((ni-i)*yt(3) + (i-1)*yt(4)) * (Nk-k) +                   &
          ((ni-i)*yt(7) + (i-1)*yt(8)) * (k-1)  ) /((ni-1)*(nk-1))
    zf5=( ((ni-i)*zt(3) + (i-1)*zt(4)) * (Nk-k) +                   &
          ((ni-i)*zt(7) + (i-1)*zt(8)) * (k-1)  ) /((ni-1)*(nk-1))
  else
    xf5 = grid % nodes(n5) % x
    yf5 = grid % nodes(n5) % y
    zf5 = grid % nodes(n5) % z
  endif

  ! Face II
  if(grid % nodes(n2) % x == HUGE) then
    xf2=( ((nj-j)*xt(1) + (j-1)*xt(3)) * (nk-k) +                   &
          ((nj-j)*xt(5) + (j-1)*xt(7)) * (k-1)  ) /((nj-1)*(nk-1))
    yf2=( ((nj-j)*yt(1) + (j-1)*yt(3)) * (nk-k) +                   &
          ((nj-j)*yt(5) + (j-1)*yt(7)) * (k-1)  ) /((nj-1)*(nk-1))
    zf2=( ((nj-j)*zt(1) + (j-1)*zt(3)) * (nk-k) +                   &
          ((nj-j)*zt(5) + (j-1)*zt(7)) * (k-1)  ) /((nj-1)*(nk-1))
  else
    xf2 = grid % nodes(n2) % x
    yf2 = grid % nodes(n2) % y
    zf2 = grid % nodes(n2) % z
  endif

  ! Face IV
  if(grid % nodes(n4) % x == HUGE) then
    xf4=( ((nj-j)*xt(2) + (j-1)*xt(4)) * (nk-k) +                   &
          ((nj-j)*xt(6) + (j-1)*xt(8)) * (k-1)  ) /((nj-1)*(nk-1))
    yf4=( ((nj-j)*yt(2) + (j-1)*yt(4)) * (nk-k) +                   &
          ((nj-j)*yt(6) + (j-1)*yt(8)) * (k-1)  ) /((nj-1)*(nk-1))
    zf4=( ((nj-j)*zt(2) + (j-1)*zt(4)) * (nk-k) +                   &
          ((nj-j)*zt(6) + (j-1)*zt(8)) * (k-1)  ) /((nj-1)*(nk-1))
  else
    xf4 = grid % nodes(n4) % x
    yf4 = grid % nodes(n4) % y
    zf4 = grid % nodes(n4) % z
  endif

  if( grid % nodes(n) % x == HUGE ) then
    grid % nodes(n) % x = ( xf1*(nk-k) + xf6*(k-1) ) * wx16 / (nk-1)        &
                 + ( xf2*(ni-i) + xf4*(i-1) ) * wx24 / (ni-1)               &
                 + ( xf3*(nj-j) + xf5*(j-1) ) * wx35 / (nj-1) 

    grid % nodes(n) % y = ( yf1*(nk-k) + yf6*(k-1) ) * wy16 / (nk-1)        &
                 + ( yf2*(ni-i) + yf4*(i-1) ) * wy24 / (ni-1)               &
                 + ( yf3*(nj-j) + yf5*(j-1) ) * wy35 / (nj-1) 

    grid % nodes(n) % z = ( zf1*(nk-k) + zf6*(k-1) ) * wz16 / (nk-1)        &
                 + ( zf2*(ni-i) + zf4*(i-1) ) * wz24 / (ni-1)               &
                 + ( zf3*(nj-j) + zf5*(j-1) ) * wz35 / (nj-1) 
  end if

  end subroutine Laplac
