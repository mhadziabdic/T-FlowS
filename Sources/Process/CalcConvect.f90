!======================================================================!
  subroutine CalcConvect
!----------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c1, c2, s
!======================================================================!

  call CalcFlux   

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- on the boundary perform the extrapolation
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
    call GraPhi(T % n,1,PHIx,.TRUE.)     ! dT/dx
    call GraPhi(T % n,2,PHIy,.TRUE.)     ! dT/dy
    call GraPhi(T % n,3,PHIz,.TRUE.)     ! dT/dz
    call GraCorNew(T % n,PHIx,PHIy,PHIz) ! needed ?
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

!---- on the boundary perform the extrapolation
      if(c2  < 0) then
        if( (TypeBC(c2) == CONVECT) ) then
          T % n(c2) = T % n(c2) - ( Ubulk(material(c1)) * PHIx(c1)        & 
                                  + Vbulk(material(c1)) * PHIy(c1)        &
                                  + Wbulk(material(c1)) * PHIz(c1) ) * dt
        end if
      end if
    end do
  end if

  RETURN 

  end subroutine CalcConvect
