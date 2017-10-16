!======================================================================!
  subroutine UserPlainGen_Tecplot_AHill(n)
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
!-----------------------------[Arguments]------------------------------!
    integer             :: n
    real                :: x1_dis, x0_dis, y1_dis, y0_dis, z1_dis, z0_dis
    real                :: TauWup, TauWdown
    real    :: Utan, UnorSq, Unor, UtotSq, dely, Stot, R, UtanX
!-------------------------------[Locals]-------------------------------!
    integer             :: c, s, c1, c2
    character*24        :: inflowfile
    character*80        :: namout
    character*39        :: path
!======================================================================!
   TauWdown = 0.0
   TauWup   = 0.0 
  
!  open(19, file='Cyl_data.dat', status='OLD', position = 'append')

  open(9, file='Slice.dat')
  if(this < 2) write(*,*) '# Now reading the file: Slice.dat '
  read(9,*) x1_dis, x0_dis
  read(9,*) y1_dis, y0_dis
  read(9,*) z1_dis, z0_dis
  if(this < 2) write(*,*) '# X:[ ', x0_dis, " ;", x1_dis, "]", &
  ' Y:[ ', y0_dis, " ;", y1_dis, "]", ' Z:[ ', z0_dis, " ;", z1_dis, "]"

    inflowfile = 'TecPla-xxxxxx-xxx.dat'

    write(inflowfile(8:13),'(I6.6)') n
    write(inflowfile(15:17),'(I3.3)') this
    open(500+this, file=inflowfile)

    if( this < 2 ) then
      write(*,*) 'Capturing field..'
    end if


!     call GradP(P % n,Px,Py,Pz)
!     call CalcShear(U%n,V%n,W%n,Shear)
!     call CalcVort()


!    call GraPhi(U % mean, 1, Ux,.TRUE.)    ! dU/dx
!    call GraPhi(U % mean, 2, Uy,.TRUE.)    ! dU/dy
!    call GraPhi(U % mean, 3, Uz,.TRUE.)    ! dU/dz

!    call GraPhi(V % mean, 1, Vx,.TRUE.)    ! dV/dx
!    call GraPhi(V % mean, 2, Vy,.TRUE.)    ! dV/dy
!    call GraPhi(V % mean, 3, Vz,.TRUE.)    ! dV/dz

!    call GraPhi(W % mean, 1, Wx,.TRUE.)    ! dW/dx
!    call GraPhi(W % mean, 2, Wy,.TRUE.)    ! dW/dy
!    call GraPhi(W % mean, 3, Wz,.TRUE.)    ! dW/dz

    do c=1, NC
      if( xc(c) < x1_dis .and. xc(c) > x0_dis .and. yc(c) < y1_dis .and. &
!        yc(c) > y0_dis .and. zc(c) < z1_dis .and. zc(c) > z0_dis) then
        yc(c) > y0_dis .and. WallDs(c) < z1_dis .and. WallDs(c) > z0_dis) then
        if (SIMULA == ZETA) then
           write (500+this,'(7E17.7E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), Kin % n(c), zc(c)
!           write (500+this,'(19E17.7E3)') xc(c), yc(c), U % mean(c), V % mean(c), W % mean(c), P % mean(c), &
!                                                  UU % mean(c) - U % mean(c)*U % mean(c), & !uu resolved
!                                                  VV % mean(c) - V % mean(c)*V % mean(c), & !vv resolved
!                                                  Ww % mean(c) - W % mean(c)*W % mean(c), & !ww resolved
!                                                  UV % mean(c) - U % mean(c)*V % mean(c), & !uv resolved
!                                                  UW % mean(c) - U % mean(c)*W % mean(c), & !uw resolved
!                                                  VW % mean(c) - V % mean(c)*W % mean(c), & !vw resolved
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Ux(c) ), & !uu modelled
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Vy(c) ), & !vv modelled
!                                                  (2.0/3.0)* Kin % mean(c) - 2.0*VISt_mean(c)*( Wz(c) ), & !ww modelled
!                                                  - VISt_mean(c)*( Uy(c) + Vx(c) ), & !uv modelled
!                                                  - VISt_mean(c)*( Uz(c) + Wx(c) ), & !uw modelled
!                                                  - VISt_mean(c)*( Vz(c) + Wy(c) ), & !vw modelled
!                                                  Eps % mean(c)
        elseif(SIMULA == K_EPS) then
           write (500+this,'(10E17.7E3)') xc(c), yc(c), zc(c), U%n(c), V%n(c), W%n(c), P%n(c), Kin%n(c),&
                                          Eps%n(c), sqrt(Kin%n(c))*Cmu25*WallDs(c)/VISc  
        elseif (SIMULA == LES) then
           write (500+this,'(6E17.7E3)') xc(c), zc(c), U % n(c), V % n(c), W%n(c), P % n(c)!, &
!           write (500+this,'(11E17.7E3)') xc(c), yc(c), U % n(c), V % n(c), W%n(c), P % n(c), &
!                                        (uu % mean(c) - U%mean(c)*U%mean(c)),&
!                                        (vv % mean(c) - V%mean(c)*V%mean(c)),&
!                                        (ww % mean(c) - W%mean(c)*W%mean(c)),&
!                                        (uv % mean(c) - U%mean(c)*V%mean(c)),&
!                                        0.5*(Shear(c)-Vort(c))
        end if
      endif
    enddo



!    do s=1,NS
!      c1=SideC(1,s)
!      c2=SideC(2,s)

!      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
!        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!
!          UtotSq = U % n(c1) * U % n(c1) &
!                 + V % n(c1) * V % n(c1) &
!                 + W % n(c1) * W % n(c1)
!          Unor = ( U % n(c1) * Sx(s)     &
!                 + V % n(c1) * Sy(s)    &
!                 + W % n(c1) * Sz(s) )   &
!               / sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))
!          UnorSq = Unor*Unor
!          if( UtotSq   >  UnorSq) then
!            Utan = sqrt(UtotSq - UnorSq)
!          else
!            Utan = TINY
!          end if

!          if(yc(c2)>0.0.and.xc(c2)<0.0) then
!            TauWup = TauWup + VISc*Utan/WallDs(c1)*sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + &
!            Sz(s)*Sz(s))
!          else if(yc(c2)<0.0.and.xc(c2)<0.0) then
!            TauWdown = TauWdown + VISc*Utan/WallDs(c1)*sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + &
!            Sz(s)*Sz(s))
!          end if         
!        end if
!      end if
!    end do  

!    call GloSum(TauWup)
!    call GloSum(TauWdown)
    call wait
     
    close(9)
    close(500+this)
    call wait
!    if ( this < 2 ) then
!      write(*,*) 'It is done'
!      write(19,*) n, TauWdown, TauWup 
!    end if

  end subroutine UserPlainGen_Tecplot_AHill
