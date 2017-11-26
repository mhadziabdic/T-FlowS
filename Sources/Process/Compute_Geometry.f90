!==============================================================================!
  subroutine Compute_Geometry(grid)
!------------------------------------------------------------------------------!
!   Calculates additional geometrical quantities.                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, m
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, ax_t, ay_t, az_t
!==============================================================================!

  !----------------------------------------------!
  !   Calculate total surface of the cell face   ! 
  !----------------------------------------------!
  do s = 1, grid % n_faces
    Scoef(s) = (  grid % sx(s)*grid % sx(s)  &
                + grid % sy(s)*grid % sy(s)  &
                + grid % sz(s)*grid % sz(s) )  
  end do

  !-------------------------------------------------------!
  !   Calculate the distance between neighbouring cells   !
  !-------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    xc1=grid % xc(c1)
    yc1=grid % yc(c1)
    zc1=grid % zc(c1) 

    xc2=grid % xc(c2) + grid % dx(s)
    yc2=grid % yc(c2) + grid % dy(s)
    zc2=grid % zc(c2) + grid % dz(s)

    grid % dx(s) = xc2-xc1
    grid % dy(s) = yc2-yc1
    grid % dz(s) = zc2-zc1
    Scoef(s) =    Scoef(s)     &
             / (  grid % dx(s)*grid % sx(s)  &
                + grid % dy(s)*grid % sy(s)  &
                + grid % dz(s)*grid % sz(s)) 
  end do  ! sides 

  !----------------------------------------------------------!
  !   Calculate interpolation coefficients for fluid phase   !
  !----------------------------------------------------------!
  do s = 1, grid % n_faces                        ! 2mat
    c1 = grid % faces_c(1,s)                                 ! 2mat
    c2 = grid % faces_c(2,s)                                 ! 2mat
                                                  ! 2mat
    fF(s) = grid % f(s)                           ! 2mat
                                                  ! 2mat 
    if( StateMat(material(c1))==FLUID .and.  &    ! 2mat
        StateMat(material(c2))==SOLID ) then      ! 2mat
      fF(s) = 1.0                                 ! 2mat
    end if                                        ! 2mat
                                                  ! 2mat
    if( StateMat(material(c1))==SOLID .and.  &    ! 2mat 
        StateMat(material(c2))==FLUID ) then      ! 2mat 
      fF(s) = 0.0                                 ! 2mat
    end if                                        ! 2mat
                                                  ! 2mat
    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then    ! 2mat
      fF(s) = 1.0                                 ! 2mat
    end if                                        ! 2mat
  end do                                          ! 2mat

  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  IsNearWall = .FALSE.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
        IsNearWall(c1) = .TRUE.  
      end if
    end if 

    if( StateMat(material(c1))==FLUID .and.  &    ! 2mat
        StateMat(material(c2))==SOLID ) then      ! 2mat
      IsNearWall(c1) = .TRUE.                     ! 2mat
    end if                                        ! 2mat
                                                  ! 2mat
    if( StateMat(material(c1))==SOLID .and.  &    ! 2mat 
        StateMat(material(c2))==FLUID ) then      ! 2mat 
      IsNearWall(c2) = .TRUE.                     ! 2mat
    end if                                        ! 2mat
  end do  ! sides 


  !-----------------------------------------------------!
  !   Calculate total surface of the monitoring plane   !
  !-----------------------------------------------------!
  do m = 1, grid % n_materials

    bulk(m) % area_x = 0.0
    bulk(m) % area_y = 0.0
    bulk(m) % area_z = 0.0

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  > 0 .or. c2 < 0.and.TypeBC(c2) == BUFFER) then
        if( (material(c1) == m) .and. (material(c1) == material(c2)) ) then
          xc1=grid % xc(c1)
          yc1=grid % yc(c1)
          zc1=grid % zc(c1)
          xc2=grid % xc(c1) + grid % dx(s)
          yc2=grid % yc(c1) + grid % dy(s)
          zc2=grid % zc(c1) + grid % dz(s)
 
          ax_t = abs(grid % sx(s))
          ay_t = abs(grid % sy(s))
          az_t = abs(grid % sz(s))

          !-------!
          !   x   !
          !-------!
          if( (xc1 <= xp(m)).and.(xc2 > xp(m)) .or. & 
              (xc2 <= xp(m)).and.(xc1 > xp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              bulk(m) % area_x = bulk(m) % area_x + 0.5 * ax_t    
            else
              bulk(m) % area_x = bulk(m) % area_x + ax_t    
            end if
          end if

          !-------!
          !   y   !
          !-------!
          if( (yc1 <= yp(m)).and.(yc2 > yp(m)) .or. & 
              (yc2 <= yp(m)).and.(yc1 > yp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              bulk(m) % area_y = bulk(m) % area_y + 0.5 * ay_t    
            else
              bulk(m) % area_y = bulk(m) % area_y + ay_t    
            end if
          end if

          !-------!
          !   z   !
          !-------!
          if( (zc1 <= zp(m)).and.(zc2 > zp(m)) .or. & 
              (zc2 <= zp(m)).and.(zc1 > zp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              bulk(m) % area_z = bulk(m) % area_z + 0.5 * az_t    
            else
              bulk(m) % area_z = bulk(m) % area_z + az_t    
            end if
          end if

        end if  ! material 1&2
      end if  ! (c2 > 0 .and. ... )
    end do

    call glosum(bulk(m) % area_x)
    call glosum(bulk(m) % area_y)
    call glosum(bulk(m) % area_z)

  end do ! m

  end subroutine
