!==============================================================================!
  subroutine TopSys(rrun) 
!------------------------------------------------------------------------------!
!   Determines the topology of the cells, faces and boundary cells.            !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: rrun
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, i
  integer :: c1, c2, m, pass
  integer :: lfn(6,4)
!==============================================================================!

  data    lfn / 1, 1, 2, 4, 3, 5,  &
                2, 5, 6, 8, 7, 7,  &
                4, 6, 8, 7, 5, 8,  &
                3, 2, 4, 3, 1, 6  /

  write(*,*) '# Now determining the topology. This may take a while !'

  !------------------------------!
  !                              !
  !   Count the boundary cells   !
  !                              !
  !------------------------------!
  NbC = 0
  do c=1,NC
    do m=1,24   ! neighbour cells
      if(grid % cells_c(m,c)  < 0) then
        NbC   = NbC + 1

        ! Remember the boundary marker, take positive value for marker
        BCmark(-NbC) =  -grid % cells_c(m,c)  

        ! Put new boundary cell into place  
        grid % cells_c(m,c)  = -NbC

        ! Material marker
        material(-NbC) = material(c)

      end if 
    end do
  end do 

  !---------------------------!
  !   Create the array with   ! 
  !   information on  sides   !
  !---------------------------!
  NS = 0     ! initialize the number of sides
  do pass = 1,3
    if(pass == 2) WallFacFst = NS+1
  do c1=1,NC
    do m=1,24 ! through all the neighbouring cells
      c2=grid % cells_c(m,c1)
      if( (pass==1).and.(c2>c1).and.(material(c1)==material(c2)) .or. &
          (pass==2).and.(c2>c1).and.(material(c1)/=material(c2)) .or. &
          (pass==3).and.(c2<0) ) then
        NS=NS+1

        ! Which volumes are connected with side NS
        SideC(1,NS)=c1
        SideC(2,NS)=c2 

        ! Which is c2 neighbour of c1 and vice versa
        do c=1,24
          if(grid % cells_c(c,c1) == c2) then 
            SideCc(NS,1)=c
          end if

          if(c2 > 0) then
            if(grid % cells_c(c,c2) == c1) then 
              SideCc(NS,2)=c
            end if 
          end if
        end do

        ! Nodes of a side NS
        if(c2  > 0) then
          if(level(c2)  > level(c1)) then
            grid % faces_n(1,NS) = grid % cells_n( lfn(SideCc(NS,2),4), c2 )
            grid % faces_n(2,NS) = grid % cells_n( lfn(SideCc(NS,2),3), c2 )
            grid % faces_n(3,NS) = grid % cells_n( lfn(SideCc(NS,2),2), c2 )
            grid % faces_n(4,NS) = grid % cells_n( lfn(SideCc(NS,2),1), c2 )
          else
            grid % faces_n(1,NS) = grid % cells_n( lfn(m,1), c1 )
            grid % faces_n(2,NS) = grid % cells_n( lfn(m,2), c1 )
            grid % faces_n(3,NS) = grid % cells_n( lfn(m,3), c1 )
            grid % faces_n(4,NS) = grid % cells_n( lfn(m,4), c1 )
          end if
        else
          grid % faces_n(1,NS) = grid % cells_n( lfn(m,1), c1 )
          grid % faces_n(2,NS) = grid % cells_n( lfn(m,2), c1 )
          grid % faces_n(3,NS) = grid % cells_n( lfn(m,3), c1 )
          grid % faces_n(4,NS) = grid % cells_n( lfn(m,4), c1 )
        end if 

      end if
    end do   ! m
  end do     ! c1
  end do     ! pass
  WallFacLst = NS

  write(*,*) '# Wall and interface faces start at: ', WallFacFst 
  write(*,*) '# Wall and interface faces end at  : ', WallFacLst

  if(.NOT. rrun) then
  NbC = 0
  do c=1,NC
    do m=1,24   ! neighbour cells
      if(grid % cells_c(m,c)  < 0) then
        NbC   = NbC + 1

        ! Restore the boundary marker, take positive value for marker
        grid % cells_c(m,c)  = -BCmark(-NbC)
      end if 
    end do
  end do 
  end if

  if(rrun) then

  !-------------------------------------------!
  !   Find the side oposite on the boundary   !
  !-------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2  < 0) then
      do i=1,6
        if(grid % cells_c(i,c1)  ==  c2) then
          if(i == 1) SideC(0,s) = grid % cells_c(6,c1)
          if(i == 2) SideC(0,s) = grid % cells_c(4,c1)
          if(i == 3) SideC(0,s) = grid % cells_c(5,c1)
          if(i == 4) SideC(0,s) = grid % cells_c(2,c1)
          if(i == 5) SideC(0,s) = grid % cells_c(3,c1)
          if(i == 6) SideC(0,s) = grid % cells_c(1,c1)
        end if
      end do
    end if
  end do 

  !-------------------------------------------------------!
  !   Find the boundary cell on which you should copy     ! 
  !-------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0 .and. CopyC(c1) /= 0) then
      if(BCmark(c2) == copy_cond(1,0)) then
        CopyC(c2) = CopyC(c1)
      end if
    end if
  end do

  end if ! rrun

  !--------------------------------------!
  !   Is there enough allocated memory   !
  !--------------------------------------!
  if( NS  > grid % max_n_faces ) then
    write(*,*) '# Error message from Generator'
    write(*,*) '# The number sides is: ', NS
    write(*,*) '# There is space available only for:', grid % max_n_faces
    write(*,*) '# Increase the number of faces in .dom file'
    write(*,*) '# and recompile the code. Good Luck !'
    stop
  end if 

  end subroutine TopSys
