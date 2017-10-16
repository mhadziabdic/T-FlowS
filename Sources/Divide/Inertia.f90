!======================================================================!
  subroutine Inertia(sub) 
!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer :: sub                           ! subdomain 
!-------------------------------[Locals]-------------------------------!
  integer :: i, n_cells_sub
  real    :: xm, ym, zm
  real    :: i_matrix(3,3), d(3), v(3,3), d_max(3)    
!======================================================================!

  xm=0.0
  ym=0.0
  zm=0.0
  n_cells_sub=0
  do i=1,NC
    if(proces(i)==sub) then
      xm = xm + xc(i)
      ym = ym + yc(i)
      zm = zm + zc(i)
      n_cells_sub=n_cells_sub+1 
    end if
  end do 
  xm=xm/n_cells_sub
  ym=ym/n_cells_sub
  zm=zm/n_cells_sub

  write(*,*) 'Center of mass for subdomain ', sub, ' is: ', xm, ym, zm

  i_matrix = 0.
  do i=1,NC
    if(proces(i)==sub) then
      i_matrix(1,1) = i_matrix(1,1) + (yc(i)-ym)**2 + (zc(i)-zm)**2
      i_matrix(2,2) = i_matrix(2,2) + (xc(i)-xm)**2 + (zc(i)-zm)**2
      i_matrix(3,3) = i_matrix(3,3) + (xc(i)-xm)**2 + (yc(i)-ym)**2

      i_matrix(1,2) = i_matrix(1,2) - (xc(i)-xm)*(yc(i)-ym) 
      i_matrix(1,3) = i_matrix(1,3) - (xc(i)-xm)*(zc(i)-zm) 
      i_matrix(2,3) = i_matrix(2,3) - (yc(i)-ym)*(zc(i)-zm) 
    end if
  end do 
  i_matrix(2,1) = i_matrix(1,2)
  i_matrix(3,1) = i_matrix(1,3)
  i_matrix(3,2) = i_matrix(2,3)

  call Compute_Eigenvalues(i_matrix, 3, 3, d, v, i)

  write(*,*) 'd=',d(1), d(2), d(3)

  write(*,*) 'v=', (v(1,i), i=1,3)
  write(*,*) '  ', (v(2,i), i=1,3)
  write(*,*) '  ', (v(3,i), i=1,3)

  if(min(d(1),d(2),d(3)) == d(1)) then
    d_max(1)=v(1,1)
    d_max(2)=v(2,1)
    d_max(3)=v(3,1)
  else if(min(d(1),d(2),d(3)) == d(2)) then
    d_max(1)=v(1,2)
    d_max(2)=v(2,2)
    d_max(3)=v(3,2)
  else if(min(d(1),d(2),d(3)) == d(3)) then
    d_max(1)=v(1,3)
    d_max(2)=v(2,3)
    d_max(3)=v(3,3)
  end if 

  write(*,*) '# Sorting the cells'
  do i=1,NC
    iin(i) = i
    criter(i) = xc(i)*d_max(1) + yc(i)*d_max(2) + zc(i)*d_max(3)
  end do
  call RISort(criter(1),iin(1),NC,2)

  return

  end subroutine Inertia
