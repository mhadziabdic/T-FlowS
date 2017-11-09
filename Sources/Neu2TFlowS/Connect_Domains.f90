!==============================================================================!
  subroutine Connect_Domains(grid)
!------------------------------------------------------------------------------!
!   Connects two problem domains, one with periodic streamwise boundary        !
!   conditions and another with inflow-outflow.                                !
!                                                                              !
!   Note:                                                                      !
!                                                                              !
!   Situations like the on depicted bellow are now working.                    !
!                                                                              !
!   +-------+   +_                                                             !
!   |       |   | ~-_                                                          !
!   |   o---|   |    ~-_                                                       !
!   |       |   |---o   +                                                      !
!   +-------+   +_      |                                                      !
!                 ~-_   |                                                      !
!                    ~-_|                                                      !
!                       +                                                      !
!                                                                              !
!   Constrains:                                                                !
!                                                                              !
!   First domain must be the channel-like, periodic in streamwise direction.   !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, c1, c11, c12, c21, c22, s1, s2
  integer           :: copy_marker,  x_copy, y_copy, z_copy
  real              :: xc_12, xc_22
  real              :: yc_12, yc_22
  real              :: zc_12, zc_22
  real,parameter    :: one3rd = 0.3333333
  real,parameter    :: two3rd = 0.6666666
  character(len=80) :: answer
!==============================================================================!

  n_copy = 0
  x_copy = 0
  y_copy = 0
  z_copy = 0

1 write(*,*) '#======================================'
  write(*,*) '# Enter the copy marker (skip to exit):'
  write(*,*) '#--------------------------------------'
  call ReadC(5,inp,tn,ts,te)
  read(inp, *) answer
  call To_Upper_Case(answer)
  if(answer == 'SKIP') then
    return
  else 
    read(inp,*) copy_marker 
  end if    

  !-------!
  !   X   !
  !-------! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) ((one3rd*s1)/(1.0*NS)) * 100, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(grid % dx(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == copy_marker) then

            yc_12 = 0.0
            zc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              yc_12 = yc_12 + grid % yn(grid % faces_n(i,s1))
              zc_12 = zc_12 + grid % zn(grid % faces_n(i,s1))
            end do
            yc_12 = yc_12 / (real(grid % faces_n_nodes(s1)))
            zc_12 = zc_12 / (real(grid % faces_n_nodes(s1)))

            yc_22 = 0.0
            zc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              yc_22 = yc_22 + grid % yn(grid % faces_n(i,s2))
              zc_22 = zc_22 + grid % zn(grid % faces_n(i,s2))
            end do
            yc_22 = yc_22 / (real(grid % faces_n_nodes(s2)))
            zc_22 = zc_22 / (real(grid % faces_n_nodes(s2)))
              
            if( Approx( yc_22, yc_12, tol=1.e-4 ) .and. &
                Approx( zc_22, zc_12, tol=1.e-4 ) ) then
              n_copy = n_copy + 1
              x_copy = x_copy + 1
              if( abs(grid % xc(c11)-grid % xc(c22)) <  &
                  abs(grid % xc(c12)-grid % xc(c22))) c1 = c11
              if( abs(grid % xc(c11)-grid % xc(c22)) >  &
                  abs(grid % xc(c12)-grid % xc(c22))) c1 = c12
              CopyS(1, n_copy) = c1
              CopyS(2, n_copy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  !-------!
  !   Y   !
  !-------! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) (one3rd + (one3rd*s1)/(1.0*NS)) * 100.0, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(grid % dy(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == copy_marker) then

            xc_12 = 0.0
            zc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              xc_12 = xc_12 + grid % xn(grid % faces_n(i,s1))
              zc_12 = zc_12 + grid % zn(grid % faces_n(i,s1))
            end do
            xc_12 = xc_12 / (real(grid % faces_n_nodes(s1)))
            zc_12 = zc_12 / (real(grid % faces_n_nodes(s1)))

            xc_22 = 0.0
            zc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              xc_22 = xc_22 + grid % xn(grid % faces_n(i,s2))
              zc_22 = zc_22 + grid % zn(grid % faces_n(i,s2))
            end do
            xc_22 = xc_22 / (real(grid % faces_n_nodes(s2)))
            zc_22 = zc_22 / (real(grid % faces_n_nodes(s2)))
              
            if( Approx( xc_22, xc_12, tol=1.e-4 ) .and. &
                Approx( zc_22, zc_12, tol=1.e-4 ) ) then
              n_copy = n_copy + 1 
              y_copy = y_copy + 1
              if( abs(grid % yc(c11)-grid % yc(c22)) <  &
                  abs(grid % yc(c12)-grid % yc(c22))) c1 = c11
              if( abs(grid % yc(c11)-grid % yc(c22)) >  &
                  abs(grid % yc(c12)-grid % yc(c22))) c1 = c12
              CopyS(1, n_copy) = c1
              CopyS(2, n_copy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  !-------!
  !   Z   !
  !-------! 
  do s1=1,NS
    if(mod(s1,100000)==0) then
      write(*,*) (two3rd + (one3rd*s1)/(1.0*NS)) * 100.0, 'Complete'
    end if
    c11 = SideC(1,s1)
    c12 = SideC(2,s1)
    if( abs(grid % dz(s1)) > 1.0e-9 ) then
      do s2=1,NS 
        c21 = SideC(1,s2)
        c22 = SideC(2,s2)
        if(c22 < 0) then
          if(BCmark(c22) == copy_marker) then

            yc_12 = 0.0
            xc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              yc_12 = yc_12 + grid % yn(grid % faces_n(i,s1))
              xc_12 = xc_12 + grid % xn(grid % faces_n(i,s1))
            end do
            yc_12 = yc_12 / (real(grid % faces_n_nodes(s1)))
            xc_12 = xc_12 / (real(grid % faces_n_nodes(s1)))

            yc_22 = 0.0
            xc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              yc_22 = yc_22 + grid % yn(grid % faces_n(i,s2))
              xc_22 = xc_22 + grid % xn(grid % faces_n(i,s2))
            end do
            yc_22 = yc_22 / (real(grid % faces_n_nodes(s2)))
            xc_22 = xc_22 / (real(grid % faces_n_nodes(s2)))
              
            if( Approx( yc_22, yc_12, tol=1.e-4 ) .and. &
                Approx( xc_22, xc_12, tol=1.e-4 ) ) then
              n_copy = n_copy + 1 
              z_copy = z_copy + 1
              if( abs(grid % zc(c11)-grid % zc(c22)) <  &
                  abs(grid % zc(c12)-grid % zc(c22))) c1 = c11
              if( abs(grid % zc(c11)-grid % zc(c22)) >  &
                  abs(grid % zc(c12)-grid % zc(c22))) c1 = c12
              CopyS(1, n_copy) = c1
              CopyS(2, n_copy) = c21           !   inside the domain
              CopyC(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  write(*,*) '# n copy = ', n_copy
  write(*,*) '# x copy = ', x_copy
  write(*,*) '# x copy = ', y_copy
  write(*,*) '# x copy = ', z_copy

  goto 1

  end subroutine
