!==============================================================================!
  subroutine Find_Sides
!------------------------------------------------------------------------------!
!  Creates the "SideC structure"                                               !
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use neu_mod 
  use gen_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
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

  very_big = max(NN,NC)

  allocate(face_coor(NC*6)); face_coor = NN*1E+30
  allocate(face_cell(NC*6)); face_cell = 0    
  allocate(Starts(NC*6)); Starts = 0    
  allocate(Ends(NC*6));   Ends = 0    
! allocate(CellC(NC,6));  CellC = 0

  !---------------------------------------------------!
  !   Fill the generic coordinates with some values   !
  !---------------------------------------------------!
  do c=1,NC
    if(grid % cells(c) % n_nodes == 4) fn = f4n
    if(grid % cells(c) % n_nodes == 5) fn = f5n
    if(grid % cells(c) % n_nodes == 6) fn = f6n
    if(grid % cells(c) % n_nodes == 8) fn = f8n 
    do j=1,6
      if(BCtype(c,j) == 0) then 

        n_f_nod = 0
        f_nod = -1
        do k=1,4
          if(fn(j,k) > 0) then
            f_nod(k) = grid % cells(c) % n( fn(j,k))
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
  call Sort_Real_By_Index(face_coor,face_cell,NC*6,2)

  !------------------------------------------------!
  !   Anotate cell faces with same coordinates     !
  !   (I am afraid that this might be influenced   !
  !      by the numerical round-off errors)        !
  !------------------------------------------------!
  Nuber = 1
  Starts(1) = 1
  do c=2,NC*6
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
  do n3=1,Nuber
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
            do n1=1,grid % cells(c1) % n_nodes
              do n2=1,grid % cells(c2) % n_nodes
                if(grid % cells(c1) % n(n1)==grid % cells(c2) % n(n2)) then
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
              if(grid % cells(c1) % n_nodes == 4) fn = f4n
              if(grid % cells(c1) % n_nodes == 5) fn = f5n
              if(grid % cells(c1) % n_nodes == 6) fn = f6n
              if(grid % cells(c1) % n_nodes == 8) fn = f8n
              do j=1,6
                if(   grid % cells(c1) % c(j) == 0  .and.   & ! not set yet
                    ( max( MatchNodes(fn(j,1)),0 ) + &
                      max( MatchNodes(fn(j,2)),0 ) + &
                      max( MatchNodes(fn(j,3)),0 ) + &
                      max( MatchNodes(fn(j,4)),0 ) == Nmatch ) ) then
                  NS = NS + 1 
                  SideC(1,NS) = c1
                  SideC(2,NS) = c2
                  SideN(NS,0) = Nmatch 
                  do k=1,4
                    if(fn(j,k) > 0) then
                      SideN(NS,k) = grid % cells(c1) % n(fn(j,k))                
                    end if
                  end do
                  grid % cells(c1) % c(j) = 1 !  -> means: set
                end if
              end do
            end if   ! Nmatch /= 2
          end if   ! c1 /= c2
        end do   ! i2
      end do   ! i1
    end if
  end do    ! do n3

  write(*,*) '# Find_Sides: Number of sides: ', NS, NS

  end subroutine Find_Sides
