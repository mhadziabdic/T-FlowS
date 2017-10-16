!==============================================================================!
  subroutine Load_Cns
!------------------------------------------------------------------------------!
! Reads:  name.cns
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, s
  character(len=80) :: dum_s, name_in
!==============================================================================!

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.cns'

  open(9, file=name_in,FORM='unformatted')
  write(*,'(A24,A)') '# Now reading the file: ', name_in

  ! Number of cells
  read(9) NC
  read(9) NbC
  read(9) NS
  read(9) NSsh
  read(9) Nmat

  allocate (material(-NbC:NC))
  read(9) (material(c), c=1,NC)
  read(9) (material(c), c=-1,-NBC,-1)

  ! Faces
  allocate (SideC(0:2,NS+NSsh))
  read(9) (SideC(0,s), s=1,NS)
  read(9) (SideC(1,s), s=1,NS)
  read(9) (SideC(2,s), s=1,NS)

  ! Boundary cells
  allocate (bcmark(-NbC-1:-1)); bcmark=0
  allocate (CopyC(-NbC:-1));    CopyC=0
  read(9) (bcmark(c), c=-1,-NbC, -1)
  read(9) (CopyC(c), c=-1,-NbC, -1)
 
  read(9) n_copy
  write(*,*) n_copy
  allocate (CopyS(2,n_copy))
  read(9) (CopyS(1,s), s=1,n_copy)
  read(9) (CopyS(2,s), s=1,n_copy)

  close(9)

  end subroutine Load_Cns
