!==============================================================================!
  subroutine Allocate_Additional(grid)
!------------------------------------------------------------------------------!
!   Allocates additional memory for Divisor                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n
  character(len=80) :: dum_s, name_in
!==============================================================================!

  ! Variables declared in gen_mod.h90:
  allocate (NewC(-grid % n_bnd_cells-1:grid % n_cells)); NewC = 0 
  allocate (NewS( grid % n_faces));                      NewS = 0
  allocate (NewN( grid % n_nodes));                      NewN = 0 

  ! Variables declared in div.h90:
  allocate (ix(-grid % n_bnd_cells:grid % n_cells));  ix=0
  allocate (iy(-grid % n_bnd_cells:grid % n_cells));  iy=0
  allocate (iz(-grid % n_bnd_cells:grid % n_cells));  iz=0
  allocate (iin(-grid % n_bnd_cells:grid % n_cells)); iin=0
  allocate (criter(grid % n_cells));   criter=0

  allocate (proces(grid % n_cells)); proces=0
  allocate (BuSeIn(grid % n_faces)); BuSeIn=0
  allocate (BuReIn(grid % n_faces)); BuReIn=0
  allocate (BufPos(grid % n_faces)); BufPos=0

  end subroutine
