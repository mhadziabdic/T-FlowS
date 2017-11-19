!======================================================================!
  subroutine UserPerturb(grid, fac, n)
!----------------------------------------------------------------------!
!   Perturbs the flow field for the channel flow.                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer :: n
  real    :: fac
!------------------------------[Calling]-------------------------------!
  real    :: Uta, UserURms, UserVRms, UserWRms
!-------------------------------[Locals]-------------------------------!
  integer :: c, seed(1)
  real    :: randn
!======================================================================!
!   See also: Uplus, URms, VRms, WRms, perturb                         !
!----------------------------------------------------------------------!

!  seed(1) = This*100 + n
!  call random_seed(PUT = seed(1:1))    ! Set user seed

  do c=1,NC

!----------------------------------!
!      add fluctuating values      !
!----------------------------------!
    call random_number(randn)
    Uta = sqrt(   Utau(material(c))**2  &
                + Vtau(material(c))**2  &
                + Wtau(material(c))**2 )

    U % n(c)  = U % n(c)   + 2.0*(0.5-randn) *                  &
             fac * Uta * UserURms(1.0-abs(zc(c)))
    call random_number(randn)
    U % o(c)  = U % o(c)  + 2.0*(0.5-randn) *                   &
             fac * Uta * UserURms(1.0-abs(zc(c)))
    call random_number(randn)
    U % oo(c) = U % oo(c) + 2.0*(0.5-randn) *                   &
             fac * Uta * UserURms(1.0-abs(zc(c)))

    call random_number(randn)
    V % n(c)  = V % n(c)   + 2.0*(0.5-randn) *                  &
             fac * Uta * UserVRms(1.0-abs(zc(c)))
    call random_number(randn)
    V % o(c)  = V % o(c)  + 2.0*(0.5-randn) *                   &
             fac * Uta * UserVRms(1.0-abs(zc(c)))
    call random_number(randn)
    V % oo(c) = V % oo(c) + 2.0*(0.5-randn) *                   &
             fac * Uta * UserVRms(1.0-abs(zc(c)))

    call random_number(randn)
    W % n(c)  = W % n(c)   + 2.0*(0.5-randn) *                  &
             fac * Uta * UserWRms(1.0-abs(zc(c)))
    call random_number(randn)
    W % o(c)  = W % o(c)  + 2.0*(0.5-randn) *                   &
             fac * Uta * UserWRms(1.0-abs(zc(c)))
    call random_number(randn)
    W % oo(c) = W % oo(c) + 2.0*(0.5-randn) *                   &
             fac * Uta * UserWRms(1.0-abs(zc(c)))

!->>> write(*,*) U(c), V(c), W(c)

    write(LinMon0(127:138),'(A6)') ' <<-- '           

  end do

  end subroutine UserPerturb
