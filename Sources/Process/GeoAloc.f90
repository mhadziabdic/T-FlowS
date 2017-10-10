!======================================================================!
  subroutine GeoAloc
!----------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                          !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use sol_mod
  use les_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!======================================================================!

  ! Variables defined in all.h90:
  allocate (xc(-NbC:NC)); xc=0.0       
  allocate (yc(-NbC:NC)); yc=0.0   
  allocate (zc(-NbC:NC)); zc=0.0  
  allocate (Sx(NS)); Sx=0.0       
  allocate (Sy(NS)); Sy=0.0       
  allocate (Sz(NS)); Sz=0.0  
  allocate (Dx(NS)); Dx=0.0       
  allocate (Dy(NS)); Dy=0.0       
  allocate (Dz(NS)); Dz=0.0  
  allocate (xsp(NS)); xsp=0.0       
  allocate (ysp(NS)); ysp=0.0       
  allocate (zsp(NS)); zsp=0.0  
  allocate (volume(-NbC:NC)); volume=0.  
  allocate (delta(-NbC:NC));  delta=0.  
  allocate (WallDs(-NbC:NC)); WallDs=0.       
  allocate (f(NS));  f=0.0  
  allocate (fF(NS)); fF=0.0  

  allocate (a1(-NbC:NC));  a1=0. 
  allocate (a2(-NbC:NC));  a2=0.  

  ! Variables defined in sol.h90:
  allocate (D % row(NC+1));     D % row=0
  allocate (D % dia(NC));       D % dia=0
  allocate (D % sav(-NbC:NC));  D % sav=0
  allocate (D % bou(-NbC:-1));  D % bou=0
  allocate (D % pos(2,NS));     D % pos=0
  allocate (p1(-NbC:NC));       p1=0
  allocate (p2(-NbC:NC));       p2=0
  allocate (q1(-NbC:NC));       q1=0
  allocate (q2(-NbC:NC));       q2=0
  allocate (r2(NC));            r2=0
  allocate (u1(NC));            u1=0
  allocate (u2(NC));            u2=0
  allocate (v1(NC));            v1=0
  allocate (v2(NC));            v2=0
  allocate (u1_plus_q1(NC));    u1_plus_q1=0

  ! Variables defined in pro_mod.h90:
  allocate (A % row(NC+1));     A % row=0
  allocate (A % dia(NC));       A % dia=0
  allocate (A % sav(-NbC:NC));  A % sav=0
  allocate (A % bou(-NbC:-1));  A % bou=0
  allocate (A % pos(2,NS));     A % pos=0
  allocate (b(NC));             b=0

  allocate (Scoef(NS)); Scoef=0.

  allocate (xp(Nmat));    xp   =0.0
  allocate (yp(Nmat));    yp   =0.0
  allocate (zp(Nmat));    zp   =0.0
  allocate (AreaX(Nmat)); AreaX=0.0
  allocate (AreaY(Nmat)); AreaY=0.0
  allocate (AreaZ(Nmat)); AreaZ=0.0

  ! Variables defined in par_mod.h90:
  allocate (BufInd(-NbC:-1)); BufInd=0

  !??????????????????????????????????!
  ! Is there enough allocated memory !
  !??????????????????????????????????!
  ! Do something !  

  end subroutine GeoAloc
