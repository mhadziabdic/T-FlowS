!======================================================================!
  SUBROUTINE UserCalcUtang_cylinder(k) 
!----------------------------------------------------------------------!
! This program reads name.1D file created by NEU or GEN and averages    
! the results in the homogeneous direction directions.            
! The results are writen in files name_res.dat and name_res_plus.dat
!----------------------------------------------------------------------!
    USE all_mod
    USE allp_mod
    USE les_mod
    USE pro_mod
    USE par_mod
    USE rans_mod
!----------------------------------------------------------------------!
    IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
    REAL :: Ufric, Wall_near 
!-------------------------------[Locals]-------------------------------!
    INTEGER             :: s, c1, c2, k
    CHARACTER           :: path*39 
    REAL                :: Fpx1,Fpy1,Fvx1,Fvy1,Fvz1, FtotalX, FtotalY 
    REAL                :: Fpx2,Fpy2,Fvx2,Fvy2,Fvz2, Mom
!--------------------------------[CVS]---------------------------------!
!  $Id: UserCalcUtang_cylinder.f90,v 1.1 2017/08/31 22:42:35 mhadziabdic Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/User/UserCalcUtang_cylinder.f90,v $  
!======================================================================!

! Fpx, Fpy, Fpz=0, rho*(Fvx, Fvy, Fvz) - pressure and viscous forces acting on a cylinder
! FtotalX, FtotalY

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
           Stot = sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s))      
           Fpx1 = Fpx1 + P%n(c1)*Sx(s) 
           Fpy1 = Fpy1 + P%n(c1)*Sy(s)
           Fvx1 = Fvx1 + (VISc)*(2.0*Ux(c1)*Sx(s)+(Uy(c1)+Vx(c1))*Sy(s)) 
           Fvy1 = Fvy1 + (VISc)*(Vx(c1)*Sx(s)+Vy(c1)*Sy(s))
           Fvz1 = Fvz1 + (VISc)*(Wx(c1)*Sx(s)+Wy(c1)*Sy(s))
           Mom  = Mom  + (Fpx1+Fvx1)*0.5*Sy(s)/Stot + (Fpy1+Fvy1)*0.5*Sx(s)/Stot
         else if(yc(c2) < 0.0) then  
           Fpx2 = Fpx2 + P%n(c1)*Sx(s)
           Fpy2 = Fpy2 + P%n(c1)*Sy(s)
           Fvx2 = Fvx2 + (VISc)*(2.0*Ux(c1)*Sx(s)+(Uy(c1)+Vx(c1))*Sy(s))
           Fvy2 = Fvy2 + (VISc)*(Vx(c1)*Sx(s)+Vy(c1)*Sy(s))
           Fvz2 = Fvz2 + (VISc)*(Wx(c1)*Sx(s)+Wy(c1)*Sy(s))
           Mom  = Mom  + (Fpx2+Fvx2)*0.5*Sy(s)/Stot + (Fpy2+Fvy2)*0.5*Sx(s)/Stot
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
   call GloSum(Mom)

   call Wait

   Omega = Omega + Mom*dtime 
 
  END SUBROUTINE UserCalcUtang_cylinder
