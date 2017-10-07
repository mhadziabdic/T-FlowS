!======================================================================!
  subroutine Split(sub, n_parts)
!----------------------------------------------------------------------!
!   Splits the domain by a geometrical multisection technique.         !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod 
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer :: sub                           ! subdomain to be splitted
  integer :: n_parts                          ! number of new partitions
!-------------------------------[Locals]-------------------------------!
  character :: dir                         ! direction for splitting
  integer   :: c, ic, j
  integer   :: n_filled
  real      :: xmax,ymax,zmax,xmin,ymin,zmin, delx,dely,delz,dxyz 
!======================================================================!
! TESTING TESTING TESTING
!--------------------------------------------! 
!     Find the smallest moment of inertia    !
!--------------------------------------------! 
  if(division_algorithm == INERTIAL) then
    call Inertia(sub)
  end if

!-------------------------------------------------! 
!     Find the largest dimension of the domain    !
!-------------------------------------------------! 
  if(division_algorithm == COORDINATE) then
    xmax=-HUGE
    ymax=-HUGE
    zmax=-HUGE
    xmin=+HUGE
    ymin=+HUGE
    zmin=+HUGE
    do c=1,NC
      if(proces(c) == sub) then
        xmax=max(xmax, xc(c))
        ymax=max(ymax, yc(c))
        zmax=max(zmax, zc(c))
        xmin=min(xmin, xc(c))
        ymin=min(ymin, yc(c))
        zmin=min(zmin, zc(c))
      end if
    end do
    delx = xmax - xmin
    dely = ymax - ymin
    delz = zmax - zmin

    dxyz = max(delx,dely,delz)
    if(delz == dxyz) dir='z'
    if(dely == dxyz) dir='y'
    if(delx == dxyz) dir='x'
  end if

  do j=1,n_parts-1

    n_filled   = 0

!----------------------------!
!     Fill the subdomain     !
!----------------------------!

    if(division_algorithm==COORDINATE) then

      if(dir == 'x') then
        do c=1,NC
          ic=ix(c)
          if(proces(ic) == sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= subNC(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if

      if(dir == 'y') then
        do c=1,NC
          ic=iy(c)
          if(proces(ic) == sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= subNC(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if

      if(dir == 'z') then
        do c=1,NC
          ic=iz(c)
          if(proces(ic) == sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= subNC(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if
    end if

    if(division_algorithm==INERTIAL) then
      do c=1,NC
        ic=iin(c)
        if(proces(ic) == sub) then
          proces(ic) = n_sub+j
          n_filled = n_filled + 1
          if(n_filled >= subNC(sub)/(n_parts-j+1)) goto 2 
        end if
      end do
    end if

!--------------------------------! 
!     Subdomain is filled up     !
!--------------------------------! 
 2  subNC(n_sub+j) = n_filled
    subNC(sub) = subNC(sub) - n_filled

  end do  ! j

  n_sub = n_sub + n_parts - 1

  end subroutine Split
