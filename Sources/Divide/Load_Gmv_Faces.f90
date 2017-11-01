!==============================================================================!
  subroutine BCelLoa
!------------------------------------------------------------------------------!
! Reads:  NAME.faces.gmv  NAME.shadow.gmv                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, n, s, dum_i
  character(len=80) :: dum_s, name_in
!==============================================================================!

  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+10) = '.faces.gmv'
  write(*,*) '# Now reading the file:', name_in
  open(9, file=name_in)

  !-----------!    
  !   Start   !
  !-----------!    
  read(9,'(A80)') dum_s 

  !------------------!    
  !   Node section   !
  !------------------!    
  read(9,'(A80)') dum_s 
  do n=1,3*NN
    read(9,'(A80)') dum_s 
  end do  

  !------------------!    
  !   Cell section   !
  !------------------!    
  read(9,'(A80)') dum_s 
  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    read(9,*) dum_s, dum_i
    if(dum_s == 'tri') then 
      grid % faces_n_nodes(s) = 3
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s)
    else if(dum_s == 'quad') then
      grid % faces_n_nodes(s) = 4
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s),  &
                grid % faces_n(4,s)  
    else
      write(*,*) 'Unsupported cell-face type:', dum_s
      write(*,*) 'Exiting'
      stop
    end if
  end do  

  close(9)

  !----------------------!    
  !   Read shadow file   !
  !----------------------!    
  call Name_File(0, name_in, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, file=name_in)
  write(6, *) 'Now reading the file:', name_in

  do s=NS+1,NS+NSsh
    read(9,*) grid % faces_n_nodes(s)
    if(grid % faces_n_nodes(s)==3) then
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s),  &
                SideC(1,s), SideC(2,s)
    else if(grid % faces_n_nodes(s)==4) then
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s),  &
                grid % faces_n(4,s),  &
                SideC(1,s), SideC(2,s)
    end if
  end do

  close(9)

  end subroutine BCelLoa
