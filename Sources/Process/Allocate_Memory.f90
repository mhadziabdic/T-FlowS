!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                                  !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use Grid_Mod
  use Matrix_Mod
  use Solvers_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Variables defined in all.h90:
  allocate (volume(-NbC:NC)); volume=0.  
  allocate (delta(-NbC:NC));  delta=0.  
  allocate (WallDs(-NbC:NC)); WallDs=0.       
  allocate (f(NS));  f=0.0  
  allocate (fF(NS)); fF=0.0  

  allocate (a1(-NbC:NC));  a1=0. 
  allocate (a2(-NbC:NC));  a2=0.  

  ! Variables defined in sol.h90:
  call Solvers_Mod_Allocate_Vectors(NbC, NC)
  call Matrix_Mod_Allocate(D, NbC, NC, NS)

  ! Variables defined in pro_mod.h90:
  call Matrix_Mod_Allocate(A, NbC, NC, NS)
  allocate (b(NC));  b=0

  allocate (Scoef(NS)); Scoef=0.

  allocate (xp(grid % n_materials));    xp   =0.0
  allocate (yp(grid % n_materials));    yp   =0.0
  allocate (zp(grid % n_materials));    zp   =0.0
  allocate (AreaX(grid % n_materials)); AreaX=0.0
  allocate (AreaY(grid % n_materials)); AreaY=0.0
  allocate (AreaZ(grid % n_materials)); AreaZ=0.0

  ! Variables defined in par_mod.h90:
  allocate (BufInd(-NbC:-1)); BufInd=0

  !??????????????????????????????????!
  ! Is there enough allocated memory !
  !??????????????????????????????????!
  ! Do something !  

  end subroutine
