!==============================================================================!
  subroutine Find_Faces(grid)
!------------------------------------------------------------------------------!
!  Creates the "SideC structure"                                               !
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use neu_mod 
  use gen_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c, c1, c2, n1, n2, n3, f_nod(4), n_f_nod
  integer             :: Nmatch, j, MatchNodes(-1:8) 
  integer             :: i1, i2, k, Nuber
  integer             :: fn(6,4)
  real   ,allocatable :: face_coor(:) 
  integer,allocatable :: face_cell(:), Starts(:), Ends(:) 
  real                :: very_big
!==============================================================================!

  very_big = max(grid % n_nodes,grid % n_cells)

  allocate(face_coor(grid % n_cells*6)); face_coor = grid % n_nodes*HUGE
  allocate(face_cell(grid % n_cells*6)); face_cell = 0    
  allocate(Starts(grid % n_cells*6)); Starts = 0    
  allocate(Ends(grid % n_cells*6));   Ends = 0    
! allocate(CellC(grid % n_cells,6));  CellC = 0

  !---------------------------------------------------!
  !   Fill the generic coordinates with some values   !
  !---------------------------------------------------!
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 4) fn = f4n
    if(grid % cells_n_nodes(c) == 5) fn = f5n
    if(grid % cells_n_nodes(c) == 6) fn = f6n
    if(grid % cells_n_nodes(c) == 8) fn = f8n 
    do j = 1, 6
      if(BCtype(c,j) == 0) then 

        n_f_nod = 0
        f_nod = -1
        do k = 1, 4
          if(fn(j,k) > 0) then
            f_nod(k) = grid % cells_n(fn(j,k), c)
            n_f_nod = n_f_nod + 1
          end if
        end do

        if( n_f_nod >  0 ) then
          if(f_nod(4) > 0) then
            face_coor((c-1)*6+j) =  very_big*(max(f_nod(1), f_nod(2), f_nod(3), f_nod(4)))   &
                                 +            min(f_nod(1), f_nod(2), f_nod(3), f_nod(4))
          else
            face_coor((c-1)*6+j) =  very_big*(max(f_nod(1), f_nod(2), f_nod(3)))   &
                                 +            min(f_nod(1), f_nod(2), f_nod(3))
           end if
          face_cell((c-1)*6+j) = c 
        end if 
      end if
    end do
  end do

  !--------------------------------------------------!
  !   Sort the cell faces according to coordinares   !
  !--------------------------------------------------!
  call Sort_Real_By_Index(face_coor,face_cell,grid % n_cells*6,2)

  !------------------------------------------------!
  !   Anotate cell faces with same coordinates     !
  !   (I am afraid that this might be influenced   !
  !      by the numerical round-off errors)        !
  !------------------------------------------------!
  Nuber = 1
  Starts(1) = 1
  do c=2,grid % n_cells*6
    if( face_coor(c) /= face_coor(c-1) ) then
      Nuber = Nuber + 1
      Starts(Nuber) = c
      Ends(Nuber-1) = c-1
    end if
  end do

  !-------------------------------------------!
  !                                           !
  !   Main loop to fill the SideC structure   !
  !                                           !
  !-------------------------------------------!
  do n3 = 1, Nuber
    if(Starts(n3) /= Ends(n3)) then
      do i1=Starts(n3),Ends(n3)
        do i2=i1+1,Ends(n3)
          c1 = min(face_cell(i1),face_cell(i2))
          c2 = max(face_cell(i1),face_cell(i2))
          if(c1 /= c2) then

            !------------------------------!
            !   Number of matching nodes   !
            !------------------------------!
            Nmatch     = 0
            MatchNodes = 0 
            do n1 = 1, grid % cells_n_nodes(c1)
              do n2 = 1, grid % cells_n_nodes(c2)
                if(grid % cells_n(n1,c1)==grid % cells_n(n2,c2)) then
                  Nmatch = Nmatch + 1 
                  MatchNodes(n1) = 1
                end if
              end do
            end do

            !-----------------------!
            !   general + general   ! 
            !     c1        c2      !
            !-----------------------!
            if(Nmatch > 2) then 
              if(grid % cells_n_nodes(c1) == 4) fn = f4n
              if(grid % cells_n_nodes(c1) == 5) fn = f5n
              if(grid % cells_n_nodes(c1) == 6) fn = f6n
              if(grid % cells_n_nodes(c1) == 8) fn = f8n
              do j = 1, 6
                if(grid % cells_c(j, c1) == 0  .and.   & ! not set yet
                    ( max( MatchNodes(fn(j,1)),0 ) + &
                      max( MatchNodes(fn(j,2)),0 ) + &
                      max( MatchNodes(fn(j,3)),0 ) + &
                      max( MatchNodes(fn(j,4)),0 ) == Nmatch ) ) then
                  grid % n_faces = grid % n_faces + 1 
                  grid % faces_c(1,grid % n_faces) = c1
                  grid % faces_c(2,grid % n_faces) = c2
                  grid % faces_n_nodes(grid % n_faces) = Nmatch 
                  do k = 1, 4
                    if(fn(j,k) > 0) then
                      grid % faces_n(k,grid % n_faces) = grid % cells_n(fn(j,k), c1)                
                    end if
                  end do
                  grid % cells_c(j, c1) = 1 !  -> means: set
                end if
              end do
            end if   ! Nmatch /= 2
          end if   ! c1 /= c2
        end do   ! i2
      end do   ! i1
    end if
  end do    ! do n3

  print *, '# Find_Faces: Number of faces: ', grid % n_faces, grid % n_faces

  end subroutine
