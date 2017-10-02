!======================================================================!
  subroutine Connect3 
!----------------------------------------------------------------------!
! Connects two problem domains, one with periodic streamwise boundary  !
! conditions and another with inflow-outflow.                          !
!                                                                      !
! Note:                                                                !
! ~~~~~                                                                !
!                                                                      !
! Situations like the on depicted bellow are now working.              !
!                                                                      !
!   +-------+   +_                                                     ! 
!   |       |   | ~-_                                                  ! 
!   |   o---|   |    ~-_                                               !
!   |       |   |---o   +                                              ! 
!   +-------+   +_      |                                              !
!                 ~-_   |                                              !
!                    ~-_|                                              !
!                       +                                              !
!                                                                      !
! Constrains:                                                          !
! ~~~~~~~~~~~                                                          !
!  1. First domain must be the channel-like, with peridicity in the    !
!     streamwise direction.                                            !
!                                                                      !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod 
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Calling]-------------------------------!
  interface
    logical function Approx(A,B,tol)
      real          :: A,B
      real,OPTIONAL :: tol
    end function Approx
  end interface
!-------------------------------[Locals]-------------------------------!
  integer      :: i, c1, c11, c12, c21, c22, s1, s2
  integer      :: CopyMarker,  x_copy, y_copy, z_copy
  real         :: xc_12, xc_22
  real         :: yc_12, yc_22
  real         :: zc_12, zc_22
  real,parameter :: one3rd = 0.3333333
  real,parameter :: two3rd = 0.6666666
  character*80 :: Answer
!======================================================================!

  Ncopy = 0
  x_copy = 0
  y_copy = 0
  z_copy = 0

1 write(*,*) 'Enter the copy marker (skip to exit):'
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) answer
  call ToUppr(answer)
  if(answer == 'SKIP') then
    return
  else 
    read(inp,*) CopyMarker 
  end if    

!---!
! x !
!---! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) ((one3rd*s1)/(1.0*NS)) * 100, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(Dx(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == CopyMarker) then

            yc_12 = 0.0
            zc_12 = 0.0
            do i=1,SideN(s1,0)
              yc_12 = yc_12 + y(SideN(s1,i))
              zc_12 = zc_12 + z(SideN(s1,i))
            end do
            yc_12 = yc_12 / (real(SideN(s1,0)))
            zc_12 = zc_12 / (real(SideN(s1,0)))

            yc_22 = 0.0
            zc_22 = 0.0
            do i=1,SideN(s2,0)
              yc_22 = yc_22 + y(SideN(s2,i))
              zc_22 = zc_22 + z(SideN(s2,i))
            end do
            yc_22 = yc_22 / (real(SideN(s2,0)))
            zc_22 = zc_22 / (real(SideN(s2,0)))
              
            if( Approx( yc_22, yc_12, tol=1.e-4 ) .and. &
                Approx( zc_22, zc_12, tol=1.e-4 ) ) then
              Ncopy = Ncopy + 1
              x_copy = x_copy + 1
              if( abs(xc(c11)-xc(c22)) < abs(xc(c12)-xc(c22))) c1 = c11
              if( abs(xc(c11)-xc(c22)) > abs(xc(c12)-xc(c22))) c1 = c12
              CopyS(1, Ncopy) = c1
              CopyS(2, Ncopy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do
!---!
! y !
!---! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) (one3rd + (one3rd*s1)/(1.0*NS)) * 100.0, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(Dy(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == CopyMarker) then

            xc_12 = 0.0
            zc_12 = 0.0
            do i=1,SideN(s1,0)
              xc_12 = xc_12 + x(SideN(s1,i))
              zc_12 = zc_12 + z(SideN(s1,i))
            end do
            xc_12 = xc_12 / (real(SideN(s1,0)))
            zc_12 = zc_12 / (real(SideN(s1,0)))

            xc_22 = 0.0
            zc_22 = 0.0
            do i=1,SideN(s2,0)
              xc_22 = xc_22 + x(SideN(s2,i))
              zc_22 = zc_22 + z(SideN(s2,i))
            end do
            xc_22 = xc_22 / (real(SideN(s2,0)))
            zc_22 = zc_22 / (real(SideN(s2,0)))
              
            if( Approx( xc_22, xc_12, tol=1.e-4 ) .and. &
                Approx( zc_22, zc_12, tol=1.e-4 ) ) then
              Ncopy = Ncopy + 1 
              y_copy = y_copy + 1
              if( abs(yc(c11)-yc(c22)) < abs(yc(c12)-yc(c22))) c1 = c11
              if( abs(yc(c11)-yc(c22)) > abs(yc(c12)-yc(c22))) c1 = c12
              CopyS(1, Ncopy) = c1
              CopyS(2, Ncopy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do
!---!
! z !
!---! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) (two3rd + (one3rd*s1)/(1.0*NS)) * 100.0, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(Dz(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == CopyMarker) then

            yc_12 = 0.0
            xc_12 = 0.0
            do i=1,SideN(s1,0)
              yc_12 = yc_12 + y(SideN(s1,i))
              xc_12 = xc_12 + x(SideN(s1,i))
            end do
            yc_12 = yc_12 / (real(SideN(s1,0)))
            xc_12 = xc_12 / (real(SideN(s1,0)))

            yc_22 = 0.0
            xc_22 = 0.0
            do i=1,SideN(s2,0)
              yc_22 = yc_22 + y(SideN(s2,i))
              xc_22 = xc_22 + x(SideN(s2,i))
            end do
            yc_22 = yc_22 / (real(SideN(s2,0)))
            xc_22 = xc_22 / (real(SideN(s2,0)))
              
            if( Approx( yc_22, yc_12, tol=1.e-4 ) .and. &
                Approx( xc_22, xc_12, tol=1.e-4 ) ) then
              Ncopy = Ncopy + 1 
              z_copy = z_copy + 1
              if( abs(zc(c11)-zc(c22)) < abs(zc(c12)-zc(c22))) c1 = c11
              if( abs(zc(c11)-zc(c22)) > abs(zc(c12)-zc(c22))) c1 = c12
              CopyS(1, Ncopy) = c1
              CopyS(2, Ncopy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  write(*,*) '# Ncopy = ', Ncopy
  write(*,*) '# X copy = ', x_copy
  write(*,*) '# Y copy = ', y_copy
  write(*,*) '# Z copy = ', z_copy

  goto 1

  end subroutine Connect3
