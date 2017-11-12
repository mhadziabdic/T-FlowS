!======================================================================!
  subroutine UserPerturb2(grid, fac, n, Dom)
!----------------------------------------------------------------------!
!   Perturbs the flow field for any flow.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use par_mod
  use les_mod
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Calling]-------------------------------!
  integer :: n, Dom
  real    :: fac
!-------------------------------[Locals]-------------------------------!
  integer :: c, seed(1), q
  real    :: randn
  integer, allocatable :: a(:) 
!======================================================================!

  call random_seed(size=q)
  allocate(a(q)) 
  a = This*100 + n
  call random_seed(PUT = a)    ! Set user seed

!----------------------------------!
!      add fluctuating values      !
!----------------------------------!
  do c=1,NC
    if(material(c) == Dom) then
    call random_number(randn)

!---- 10 % of the maximum values
    U % n(c)  = U % n(c) + .1*fac*(0.5-randn) & 
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    U % o(c)  = U % n(c)
    U % oo(c) = U % n(c)

!---- 10 % of the maximum values
    V % n(c)  = V % n(c)   + .1*fac*(0.5-randn) &
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    V % o(c)  = V % n(c)
    V % oo(c) = V % n(c) 

!---- 10 % of the maximum values
    W % n(c)  = W % n(c)   + .1*fac*(0.5-randn) &
              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
    W % o(c)  = W % n(c)
    W % oo(c) = W % n(c)
    end if
  end do


  end subroutine UserPerturb2
