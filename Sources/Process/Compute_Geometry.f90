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
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, AreaTx, AreaTy, AreaTz
!==============================================================================!

  !----------------------------------------------!
  !   Calculate total surface of the cell face   ! 
  !----------------------------------------------!
  do s=1,NS
    Scoef(s) = (Sx(s)*Sx(s)+Sy(s)*Sy(s)+Sz(s)*Sz(s))  
  end do

  !-------------------------------------------------------!
  !   Calculate the distance between neighbouring cells   !
  !-------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

    xc1=grid % xc(c1)
    yc1=grid % yc(c1)
    zc1=grid % zc(c1) 

    xc2=grid % xc(c2) + Dx(s)
    yc2=grid % yc(c2) + Dy(s)
    zc2=grid % zc(c2) + Dz(s)

    Dx(s) = xc2-xc1
    Dy(s) = yc2-yc1
    Dz(s) = zc2-zc1
    Scoef(s) = Scoef(s)/(Dx(s)*Sx(s) + Dy(s)*Sy(s) + Dz(s)*Sz(s)) 
  end do  ! sides 

  !----------------------------------------------------------!
  !   Calculate interpolation coefficients for fluid phase   !
  !----------------------------------------------------------!
  do s=1,NS                                       ! 2mat
    c1=SideC(1,s)                                 ! 2mat
    c2=SideC(2,s)                                 ! 2mat
                                                  ! 2mat
    fF(s) = f(s)                                  ! 2mat
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

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

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
  AreaX = 0.0
  AreaY = 0.0
  AreaZ = 0.0
  do m = 1, grid % n_materials
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2  > 0 .or. c2 < 0.and.TypeBC(c2) == BUFFER) then
        if( (material(c1) == m) .and. (material(c1) == material(c2)) ) then
          xc1=grid % xc(c1)
          yc1=grid % yc(c1)
          zc1=grid % zc(c1)
          xc2=grid % xc(c1)+Dx(s)
          yc2=grid % yc(c1)+Dy(s)
          zc2=grid % zc(c1)+Dz(s)
 
          AreaTx = abs(Sx(s))
          AreaTy = abs(Sy(s))
          AreaTz = abs(Sz(s))

          !-------!
          !   x   !
          !-------!
          if( (xc1 <= xp(m)).and.(xc2 > xp(m)) .or. & 
              (xc2 <= xp(m)).and.(xc1 > xp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              AreaX(m) = AreaX(m) + 0.5 * AreaTx    
            else
              AreaX(m) = AreaX(m) + AreaTx    
            end if
          end if

          !-------!
          !   y   !
          !-------!
          if( (yc1 <= yp(m)).and.(yc2 > yp(m)) .or. & 
              (yc2 <= yp(m)).and.(yc1 > yp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              AreaY(m) = AreaY(m) + 0.5 * AreaTy    
            else
              AreaY(m) = AreaY(m) + AreaTy    
            end if
          end if

          !-------!
          !   z   !
          !-------!
          if( (zc1 <= zp(m)).and.(zc2 > zp(m)) .or. & 
              (zc2 <= zp(m)).and.(zc1 > zp(m)) ) then

            ! Watch out: buffer cell faces will be counted twice
            if(c2  < 0.and.TypeBC(c2) == BUFFER) then 
              AreaZ(m) = AreaZ(m) + 0.5 * AreaTz    
            else
              AreaZ(m) = AreaZ(m) + AreaTz    
            end if
          end if

        end if  ! material 1&2
      end if  ! (c2 > 0 .and. ... )
    end do

    call glosum(AreaX(m))
    call glosum(AreaY(m))
    call glosum(AreaZ(m))

  end do ! m

  end subroutine
