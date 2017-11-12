!======================================================================!
  subroutine UserCylGen_Tecplot(grid, n)
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
!-----------------------------[Arguments]------------------------------!
    integer :: n
    real    :: x1_dis, x0_dis, y1_dis, y0_dis, z1_dis, z0_dis
    real    :: TauWup, TauWdown, r23, TKE, Cmu_mod, Pk_tmp, Cs_mod, Cs_prime
    real    :: Utan, UnorSq, Unor, UtotSq, dely, Stot, R, UtanX, TKE_mean
    real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32, a_mn_a_mn
    real    :: u11, u22, u33, u12, u13, u21, u31, u23, u32
    real    :: S11, S22, S33, S12, S13, S21, S31, S23, S32, a_lk_s_lk
!-------------------------------[Locals]-------------------------------!
    integer             :: c, s, c1, c2
    character*24        :: inflowfile
    character*80        :: namout
    character*39        :: path
!======================================================================!
!  open(9, file='Slice2.dat')
!  if(this < 2) write(*,*) '# Now reading the file: Slice.dat '
!  read(9,*) x1_dis, x0_dis
!  read(9,*) y1_dis, y0_dis
!  read(9,*) z1_dis, z0_dis
!  if(this < 2) write(*,*) '# X:[ ', x0_dis, " ;", x1_dis, "]", &
!  ' Y:[ ', y0_dis, " ;", y1_dis, "]", ' Z:[ ', z0_dis, " ;", z1_dis, "]"

    inflowfile = 'Cylind-xxxxxx-xxx.dat'

    write(inflowfile(8:13),'(I6.6)') n
    write(inflowfile(15:17),'(I3.3)') this
    open(500+this, file=inflowfile)

    if( this < 2 ) then
      write(*,*) 'Capturing field..'
    end if


    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          write (500+this,'(8E17.7E3)') xc(c1), yc(c1), zc(c1), U % n(c1), V % n(c1), &
          W%n(c1), P % n(c1), sqrt(U % n(c1)**2+V % n(c1)**2+W%n(c1)**2)
        end if
      end if
    end do  

    close(9)
    close(500+this)

  end subroutine UserCylGen_Tecplot
