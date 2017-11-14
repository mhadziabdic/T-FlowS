!==============================================================================!
  subroutine Save_Shadows(grid, sub, NNsub, NCsub)
!------------------------------------------------------------------------------!
! Writes .shadow file.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub, NCsub, NmaterBC
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c,  c1,  c2,  n, s
  character(len=80) :: name_out
!==============================================================================!

  !------------------------!
  !                        !
  !   Create shadow file   !
  !                        !
  !------------------------!
  if(sub /= 0) return

  call Name_File(sub, name_out, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', name_out

  do s = grid % n_faces+1, grid % n_faces+NSsh
    write(9,*) grid % faces_n_nodes(s) 
    if(grid % faces_n_nodes(s)==3) then
      write(9,*) grid % faces_n(1,s),  &
                 grid % faces_n(2,s),  &
                 grid % faces_n(3,s),  &
                 SideC(1,s), SideC(2,s) 
    else if(grid % faces_n_nodes(s)==4) then
      write(9,*) grid % faces_n(1,s),  &
                 grid % faces_n(2,s),  &
                 grid % faces_n(3,s),  &
                 grid % faces_n(4,s),  &
                 SideC(1,s), SideC(2,s) 
    end if
  end do  

  close(9)

  end subroutine
