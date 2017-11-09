!==============================================================================!
  subroutine Save_Cas(grid, sub, NNsub, NCsub, NSsub)
!------------------------------------------------------------------------------!
!   Writes: ".cas" file                                                        !
!                                                                              !
!   See also: number                                                           !
!   NSsub holds (has to hold) NS + NSsh                                        !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub, NCsub, NSsub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c,  c1,  c2,  n, s, Nfac, NtotFac 
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .cas file   !
  !                      !
  !----------------------!
  call Name_File(sub, name_out, '.cas', len_trim('.cas'))
  open(9, file=name_out)
  write(*,*) '# Now creating the file:', name_out

  !-----------!
  !           !
  !   Start   !
  !           !
  !-----------!
  write(9,'(A34)') '(0 "============================")'
  write(9,'(A34)') '(0 "Created by TFlowS-Generator ")'
  write(9,'(A34)') '(0 "============================")'

  !---------------!
  !               !
  !   Dimension   !
  !               !
  !---------------!
  write(9,'(A15)') '(0 "Dimension")'
  write(9,'(A5)')  '(2 3)'              

  !-----------!
  !           !
  !   Nodes   !
  !           !
  !-----------+-------------------------------------------!
  !   (10 (zone_id  first_index  last_index  type  ND))   !
  !-------------------------------------------------------!
  write(9,'(A11)') '(0 "Nodes")'

  ! Declaration section
  write(9,'(A7,Z9,Z9,A4)') '(10 (0 ', 1, NNsub, ' 1))'

  ! Regular node section
  write(9,'(A7,Z9,Z9,A4)') '(10 (7 ', 1, NNsub, ' 1)('
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(3E15.7)') grid % xn(n),  &
                                          grid % yn(n),  &
                                          grid % zn(n)
  end do
  write(9,'(A2)') '))'

  !-----------!
  !           !
  !   Faces   !
  !           !
  !-----------+-----------------------------------------------------!
  !   (13 (zone_id  first_index  last_index  type  element_type))   !
  !-----------------------------------------------------------------!

  ! Declaration section
  write(9,'(A11)') '(0 "Sides")'
  write(9,'(A7,Z9,Z9,A4)') '(13 (0 ', 1, NSsub, ' 0))'

  ! Regular side section

  !---------------------------!
  !   Faces on the boundary   !
  !---------------------------!
  NtotFac = 0
  BCmark(-NbC-1) = 20  ! set the type for periodic.
  ! It has to be 19+1, where 19 is max number of boundary. 
  ! See ReadFluentNeu.f90

  do n=1,19   ! browse through boundary condition types
    Nfac = 0
    do s=1,NSsub   ! count the faces with boundary condition "n" 
        c2 = SideC(2,s)
        if(c2 < 0) then
          if(BCmark(c2) == n) Nfac=Nfac+1
        end if 
    end do    ! faces 

    if(Nfac /= 0) then
      write(9,'(A26,I3,A3)') '(0 "Sides on the boundary ', n, ' ")'
      write(9,'(A5,Z9,Z9,Z9,A6)') '(13 (', 100+n, NtotFac+1, NtotFac+Nfac, ' 3 0)('
      do s=1,NS+NSsh
          c1 = SideC(1,s)
          c2 = SideC(2,s)
          if(c2 < 0) then  
            if(BCmark(c2) == n) then 
              if(grid % faces_n_nodes(s) == 3) then
                write(9,'(6Z9)')                 &
                  3, NewN(grid % faces_n(1,s)),  &
                     NewN(grid % faces_n(2,s)),  &
                     NewN(grid % faces_n(3,s)),  &
                     NewC(c1), 0  ! <-- set c2 to zero 
              else if(grid % faces_n_nodes(s) == 4) then
                write(9,'(7Z9)')                 &
                  4, NewN(grid % faces_n(1,s)),  &
                     NewN(grid % faces_n(2,s)),  & 
                     NewN(grid % faces_n(3,s)),  &
                     NewN(grid % faces_n(4,s)),  &
                     NewC(c1), 0  ! <-- set c2 to zero  
              end if 
            end if
          end if
      end do
      write(9,'(A2)') '))'
    end if

    ! Prepare for next boundary
    write(*,*) '# Number of faces:', Nfac, NtotFac+1, NtotFac+Nfac
    NtotFac = NtotFac+Nfac

  end do   ! n -> boundary condition types 

  ! periodic.shadow
  write(9,'(A26,I3,A3)') '(0 "Sides on the boundary ', n, ' ")'
  write(9,'(A5,Z9,Z9,Z9,A6)') '(13 (', 5, NtotFac+1, NtotFac+NSsh/2, ' 8 0)('
  do s=NS+1,NS+NSsh,2  ! =--> or maybe: do s=NS+1,NS+NSsh/2
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if(c2 < 0) then
      if(BCmark(c2) == n) then
        if(grid % faces_n_nodes(s) == 3) then
          write(9,'(6Z9)')                 &
            3, NewN(grid % faces_n(1,s)),  &
               NewN(grid % faces_n(2,s)),  &
               NewN(grid % faces_n(3,s)),  &
               NewC(c1), 0  ! <-- set c2 to zero
        else if(grid % faces_n_nodes(s) == 4) then
          write(9,'(7Z9)')                 &
            4, NewN(grid % faces_n(1,s)),  & 
               NewN(grid % faces_n(2,s)),  &
               NewN(grid % faces_n(3,s)),  &
               NewN(grid % faces_n(4,s)),  &
               NewC(c1), 0  ! <-- set c2 to zero
        end if
      end if
    end if
  end do
  write(9,'(A2)') '))' 

  ! Periodic
  write(9,'(A26,I3,A3)') '(0 "Sides on the boundary ', n, ' ")'
  write(9,'(A5,Z9,Z9,Z9,A6)') '(13 (', 6, NtotFac+1+NSsh/2, NtotFac+NSsh, ' c 0)('
  do s=NS+2,NS+NSsh,2  ! =--> or maybe: do s=NS+1,NS+NSsh/2
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == n) then
          if(grid % faces_n_nodes(s) == 3) then
            write(9,'(6Z9)')                 &
              3, NewN(grid % faces_n(1,s)),  &
                 NewN(grid % faces_n(2,s)),  &
                 NewN(grid % faces_n(3,s)),  &
                 NewC(c1), 0  ! <-- set c2 to zero
          else if(grid % faces_n_nodes(s) == 4) then
            write(9,'(7Z9)')                 &
              4, NewN(grid % faces_n(1,s)),  &
                 NewN(grid % faces_n(2,s)),  &
                 NewN(grid % faces_n(3,s)),  &
                 NewN(grid % faces_n(4,s)),  &
                 NewC(c1), 0  ! <-- set c2 to zero
          end if
        end if
      end if
  end do
  write(9,'(A2)') '))' 

  write(9,'(A7,Z9,Z9,A6)') '(18 (', NtotFac+1, NtotFac+NSsh/2, ' 6 5)('
  do s=NtotFac+1, NtotFac+NSsh/2     ! =-> ili mozda: "do s=NS+1,NS+NSsh,2" ??
    write(9,'(Z9,Z9)') s, s+NSsh/2 
  end do
  write(9,'(A2)') '))' 

  NtotFac = NtotFac + NSsh

  !-------------------------!
  !   Faces in the domain   !
  !-------------------------!

  ! First count the cell faces on the material interface
  Nfac = 0
  do s=1,NSsub
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 > 0 .and. (material(NewC(c1)) /= material(NewC(c2))) ) then
        NFac = NFac+1 
      end if
  end do
  write(*,*) '# Number of cell faces at interface: ', Nfac

  write(9,'(A33)') '(0 "Sides on material interface")'
  write(9,'(A7,Z9,Z9,A6)') '(13 (3 ', NtotFac+1, NtotFac+Nfac, ' 2 0)('
  do s=1,NSsub
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 > 0 .and. (material(NewC(c1)) /= material(NewC(c2))) ) then
        if(grid % faces_n_nodes(s) == 3) then
          write(9,'(6Z9)')                 &
            3, NewN(grid % faces_n(1,s)),  &
               NewN(grid % faces_n(2,s)),  &
               NewN(grid % faces_n(3,s)),  &
               NewC(c1), NewC(c2)
        else if(grid % faces_n_nodes(s) == 4) then
          write(9,'(7Z9)')                 &
            4, NewN(grid % faces_n(1,s)),  &
               NewN(grid % faces_n(2,s)),  & 
               NewN(grid % faces_n(3,s)),  &
               NewN(grid % faces_n(4,s)),  &
               NewC(c1), NewC(c2)
        end if 
      end if
  end do
  write(9,'(A2)') '))'

  ! Faces in the domain
  NtotFac = NtotFac+Nfac
  write(9,'(A25)') '(0 "Sides in the domain")'
  write(9,'(A7,Z9,Z9,A6)') '(13 (4 ', NtotFac+1, NSsub-NSsh/2, ' 2 0)('
  do s=1,NSsub
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 < 0) c2=0
      if(c2 > 0 .and. (material(NewC(c1)) == material(NewC(c2))) .and. &
         grid % dx(s) == 0.0 .and.  &
         grid % dy(s) == 0.0 .and.  &
         grid % dz(s) == 0.0 ) then
        if(grid % faces_n_nodes(s) == 3) then
          write(9,'(6Z9)')                 &
            3, NewN(grid % faces_n(1,s)),  &
               NewN(grid % faces_n(2,s)),  &
               NewN(grid % faces_n(3,s)),  &
               NewC(c1), NewC(c2)
        else if(grid % faces_n_nodes(s) == 4) then
          write(9,'(7Z9)')                 &
            4, NewN(grid % faces_n(1,s)),  &
               NewN(grid % faces_n(2,s)),  &
               NewN(grid % faces_n(3,s)),  &
               NewN(grid % faces_n(4,s)),  &
               NewC(c1), NewC(c2)
        end if
      end if
  end do
  write(9,'(A2)') '))'

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------+-----------------------------------------------------!
  !   (12 (zone_id  first_index  last_index  type  element_type))   !
  !-----------------------------------------------------------------!
  write(9,'(A11)') '(0 "Cells")'

  ! Declaration section
  write(9,'(A7,Z9,Z9,A4)') '(12 (0 ', 1, NCsub, ' 0))'

  ! Regular cell section
  write(9,'(A7,Z9,Z9,A6)') '(12 (1 ', 1, NCsub, ' 1 0)('
  do c=1,NC
    if(NewC(c) /= 0) then
      if(grid % cells_n_nodes(c) == 8) then       ! hexahedra   
        write(9, *) ' 4 ' 
      else if(grid % cells_n_nodes(c) == 6) then  ! prism
        write(9, *) ' 6 '
      else if(grid % cells_n_nodes(c) == 4) then  ! tetrahedra
        write(9, *) ' 2 '
      else if(grid % cells_n_nodes(c) == 5) then  ! pyramid    
        write(9, *) ' 5 '
      else
        write(*,*) 'Unsupported cell type with ', &
                    grid % cells_n_nodes(c), ' nodes.'
        write(*,*) 'Exiting'
        stop 
      end if 
    end if
  end do  
  write(9,'(A2)') '))'

  !-----------!
  !           !
  !   Zones   !
  !           !
  !-----------!
  write(9,'(A11)') '(0 "Zones")'
  write(9,'(A24)') '(45 (1 fluid fluid) () )'
  write(9,'(A22)') '(45 (2 wall wall) () )'
  write(9,'(A40)') '(45 (3 interior material-interface) () )'
  write(9,'(A38)') '(45 (4 interior default-interior) () )'
  write(9,'(A37)') '(45 (5 shadow periodic.1-shadow) () )'
  write(9,'(A32)') '(45 (6 periodic periodic.1) () )'
  write(9,'(A5,I6,A18)') '(45 (', 101 , ' wall wall-01) ())'
  write(9,'(A5,I6,A18)') '(45 (', 102 , ' wall wall-02) ())'
  write(9,'(A5,I6,A18)') '(45 (', 103 , ' wall wall-03) ())'
  write(9,'(A5,I6,A18)') '(45 (', 104 , ' wall wall-04) ())'
  write(9,'(A5,I6,A18)') '(45 (', 105 , ' wall wall-05) ())'
  write(9,'(A5,I6,A18)') '(45 (', 106 , ' wall wall-06) ())'
  write(9,'(A5,I6,A18)') '(45 (', 107 , ' wall wall-07) ())'
  write(9,'(A5,I6,A18)') '(45 (', 108 , ' wall wall-08) ())'
  write(9,'(A5,I6,A18)') '(45 (', 109 , ' wall wall-09) ())'
  write(9,'(A5,I6,A18)') '(45 (', 110 , ' wall wall-10) ())'
  write(9,'(A5,I6,A18)') '(45 (', 111 , ' wall period ) ())'

  end subroutine Save_Cas
