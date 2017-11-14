!==============================================================================!
  subroutine CalcConvect(grid)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  call Compute_Fluxes(grid)

  do s = 1, grid % n_faces
    c1=SideC(1,s)
    c2=SideC(2,s)

    ! On the boundary perform the extrapolation
    if(c2  < 0) then
      if( (TypeBC(c2) == CONVECT) ) then
        U % n(c2) = U % n(c2) - ( Ubulk(material(c1)) * Ux(c1)        & 
                                + Vbulk(material(c1)) * Uy(c1)        &
                                + Wbulk(material(c1)) * Uz(c1) ) * dt
        V % n(c2) = V % n(c2) - ( Ubulk(material(c1)) * Vx(c1)        & 
                                + Vbulk(material(c1)) * Vy(c1)        &
                                + Wbulk(material(c1)) * Vz(c1) ) * dt
        W % n(c2) = W % n(c2) - ( Ubulk(material(c1)) * Wx(c1)        & 
                                + Vbulk(material(c1)) * Wy(c1)        &
                                + Wbulk(material(c1)) * Wz(c1) ) * dt
      end if
    end if
  end do

  if(HOT==YES) then
    call GraPhi(grid, T % n, 1, phix, .TRUE.)     ! dT/dx
    call GraPhi(grid, T % n, 2, phiy, .TRUE.)     ! dT/dy
    call GraPhi(grid, T % n, 3, phiz, .TRUE.)     ! dT/dz
    call GraCorNew(T % n,phix,phiy,phiz) ! needed ?
    do s = 1, grid % n_faces
      c1=SideC(1,s)
      c2=SideC(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (TypeBC(c2) == CONVECT) ) then
          T % n(c2) = T % n(c2) - ( Ubulk(material(c1)) * phix(c1)        & 
                                  + Vbulk(material(c1)) * phiy(c1)        &
                                  + Wbulk(material(c1)) * phiz(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
