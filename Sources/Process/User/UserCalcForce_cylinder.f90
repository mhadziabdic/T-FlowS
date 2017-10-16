!======================================================================!
  subroutine UserCalcForce_cylinder(k) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous direction directions.            
! The results are writen in files name_res.dat and name_res_plus.dat
!----------------------------------------------------------------------!
    use all_mod
    use allp_mod
    use les_mod
    use pro_mod
    use par_mod
    use rans_mod
!----------------------------------------------------------------------!
    implicit none
!-----------------------------[Arguments]------------------------------!
    real :: Ufric, Wall_near 
!-------------------------------[Locals]-------------------------------!
    integer             :: s, c1, c2, k
    character           :: path*39 
    real                :: Fpx1,Fpy1,Fvx1,Fvy1,Fvz1, FtotalX, FtotalY 
    real                :: Fpx2,Fpy2,Fvx2,Fvy2,Fvz2
!======================================================================!

! Fpx, Fpy, Fpz=0, rho*(Fvx, Fvy, Fvz) - pressure and viscous forces acting on a cylinder
! FtotalX, FtotalY

!  path ='/home/IUS/mhadziabdic/Documents/T-FlowS-Sep2009/Test_Cases/RANS/CYLINDER/LES_Feb_2016'

  open(19, file='Force_cyl.dat', status='OLD', position = 'append')
!  if(this < 2) write(19,'(A50)') 'k, Fpx, Fpy, Fvx, Fvy, Fvz, FtotalX, FtotalY'
!  call wait
!    if(Ndtt == 0) then
!      open(999,file=path//'force.dat',POSITION='APPEND')
!      write(999,*) '#*******************************************************'
!      write(999,*) '#timestep:  Fpx:  Fpy:  Fvx:  Fvy:  Fvz:  FtotX:  FtotY:'
!      write(999,*) '#*******************************************************'
!      close(999)
!    if(mod(k,10)==0) then
!      open(999,file=path//'force.dat',POSITION='APPEND')


!   call GraPhi(U % mean, 1, Ux,.TRUE.)    ! dU/dx
!   call GraPhi(U % mean, 2, Uy,.TRUE.)    ! dU/dy
!   call GraPhi(U % mean, 3, Uz,.TRUE.)    ! dU/dz
!   call GraPhi(V % mean, 1, Vx,.TRUE.)    ! dV/dx
!   call GraPhi(V % mean, 2, Vy,.TRUE.)    ! dV/dy
!   call GraPhi(V % mean, 3, Vz,.TRUE.)    ! dV/dz
!   call GraPhi(W % mean, 1, Wx,.TRUE.)    ! dW/dx
!   call GraPhi(W % mean, 2, Wy,.TRUE.)    ! dW/dy 
!   call GraPhi(W % mean, 3, Wz,.TRUE.)    ! dW/dz

   Fpx1 = 0.0
   Fpy1 = 0.0
   Fvx1 = 0.0
   Fvy1 = 0.0
   Fvz1 = 0.0
   Fpx2 = 0.0
   Fpy2 = 0.0
   Fvx2 = 0.0
   Fvy2 = 0.0
   Fvz2 = 0.0
   FtotalX = 0.0
   FtotalY = 0.0

   call wait

   do s=1,NS
     c1=SideC(1,s)
     c2=SideC(2,s)
     if(c2  < 0) then
       if(bcmark(c2) == 1) then !cylinders wall
         if(yc(c2) > 0.0) then      
           Fpx1 = Fpx1 + P%n(c1)*Sx(s)
           Fpy1 = Fpy1 + P%n(c1)*Sy(s)
           Fvx1 = Fvx1 + (VISc)*(2.0*Ux(c1)*Sx(s)+(Uy(c1)+Vx(c1))*Sy(s)) 
           Fvy1 = Fvy1 + (VISc)*(Vx(c1)*Sx(s)+Vy(c1)*Sy(s))
           Fvz1 = Fvz1 + (VISc)*(Wx(c1)*Sx(s)+Wy(c1)*Sy(s))
         else if(yc(c2) < 0.0) then  
           Fpx2 = Fpx2 + P%n(c1)*Sx(s)
           Fpy2 = Fpy2 + P%n(c1)*Sy(s)
           Fvx2 = Fvx2 + (VISc)*(2.0*Ux(c1)*Sx(s)+(Uy(c1)+Vx(c1))*Sy(s))
           Fvy2 = Fvy2 + (VISc)*(Vx(c1)*Sx(s)+Vy(c1)*Sy(s))
           Fvz2 = Fvz2 + (VISc)*(Wx(c1)*Sx(s)+Wy(c1)*Sy(s))
         end if  
       end if
     end if
   end do

   FtotalX = 0.0 + Fpx1 + Fvx1 + Fpx2 + Fvx2
   FtotalY = 0.0 + Fpy1 + Fvy1 + Fpx2 + Fvx2

   call GloSum(Fpx1)
   call GloSum(Fpy1)
   call GloSum(Fvx1)
   call GloSum(Fvy1)
   call GloSum(Fvz1)
   call GloSum(Fpx2)
   call GloSum(Fpy2)
   call GloSum(Fvx2)
   call GloSum(Fvy2)
   call GloSum(Fvz2)
   call GloSum(FtotalX)
   call GloSum(FtotalY)

   call Wait
 
!   if(this<2) write(*,*) 'in CalcForce; Fpy = ', Fpy

!   if(this < 2 .and. mod(k,10)==0) then
   if(this < 2) then
     write(19,'(I9,12E15.7)') k, Fpx1, Fpy1, Fvx1, Fvy1, Fvz1, Fpx2, Fpy2, Fvx2,&
     Fvy2, Fvz2, FtotalX, FtotalY
   endif
     
!   call Wait
! Pia,b - previous values 
!!! min and max lift force
!   if(mod(k,2)==0) then
!     if ( (Fpy > 0.0) .and. (cPa > 0.0) .and. (cPb > 0.0) .and. (cPb > Fpy) .and. (cPb > cPa) ) then ! find P1 max lift force
!       cNPmax = k                    
!       ForceMAX = Fpy
!     elseif( (Fpy < 0.0) .and. (cPa < 0.0) .and. (cPb < 0.0) .and. (cPb < cPa).and. (Fpy > cPb) ) then ! find P3 min lift force
!       cNPmin = k
!       ForceMIN = Fpy 
! forecase R1, R2
! cR1 = cNP3 + FLOOR( (cNP3 - cNP2)/3.0 )
! cR2 = cNP3 + 2*FLOOR( (cNP3 - cNP2)/3.0 )
!     endif
!! going from minimum to zero
!     if( (ForceMIN*(0.6666 + 0.03)<Fpy).and.(Fpy<ForceMIN*(0.6666 - 0.03)).and.(cPb>cPa).and.(Fpy>cPb)) then
!       cNPminA = k
!     elseif( (ForceMIN*(0.3333 + 0.03)<Fpy).and.( Fpy<ForceMIN*(0.3333 - 0.03)).and.(cPb>cPa).and.(Fpy>cPb)) then
!       cNPminB = k
!     endif 
!!! zero lift force
!     if( (Fpy < 0.0) .and. (cPa > 0.0) .and. (cPa > cPb) .and. (cPb > Fpy) ) then ! find P2 zero lift force
!        cNP0a = k
!      ! forecast L1, L2
!      ! cL1 = cNP2 + FLOOR( (cNP2 - cNP1)/3.0 )
!      ! cL2 = cNP2 + 2*FLOOR( (cNP2 - cNP1)/3.0 )
!      elseif( (Fpy > 0.0) .and. (cPa < 0.0) .and. (cPa < cPb) .and. (cPb < Fpy)) then ! find P4 zero lift force
!        cNP0b = k
!      endif 
!      !! going from zero to maximum
!      if((-ForceMIN*(0.3333-0.03)<Fpy).and.(Fpy<-ForceMIN*(0.3333+0.03)).and.(cPa<cPb).and.(cPb<Fpy)) then
!      cNPmaxA = k
!      elseif((-ForceMIN*(0.6666-0.03)<Fpy).and.(Fpy<-ForceMIN*(0.6666+0.03)).and.(cPa<cPb).and.(cPb<Fpy)) then
!      cNPmaxB = k
!      endif 
! endif !mod k, 2
!      cPa = cPb
!      cPb = Fpy
!      call Wait
!  endif

   close(19)

  end subroutine UserCalcForce_cylinder
