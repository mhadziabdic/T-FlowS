!==============================================================================!
  subroutine Balance_Mass(grid)
!------------------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the mass.            ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
  use Bulk_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: m, s, c1, c2
  real    :: fac(256) 
!==============================================================================!

  !--------------------------------------!
  !   Calculate the inflow mass fluxes   !
  !--------------------------------------!
  do m = 1, grid % n_materials
    bulk(m) % mass_in = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        Flux(s) = DENc(material(c1))*( U % n(c2)*grid % sx(s) + &
                                       V % n(c2)*grid % sy(s) + &
                                       W % n(c2)*grid % sz(s) )
        if(TypeBC(c2)  ==  INFLOW) then
          if(material(c1) == m) then 
            bulk(m) % mass_in = bulk(m) % mass_in - Flux(s)
          end if
        endif
        if(TypeBC(c2)  ==  PRESSURE.and.Flux(s) < 0.0) then
          if(material(c1) == m) then
            bulk(m) % mass_in = bulk(m) % mass_in - Flux(s)
          end if
        end if
        if(TypeBC(c2)  ==  CONVECT.and.Flux(s) < 0.0) then
          if(material(c1) == m) then
            bulk(m) % mass_in = bulk(m) % mass_in - Flux(s)
          end if
        end if
      end if
    end do
    call glosum(bulk(m) % mass_in)
  end do                    

  !---------------------------------------!
  !   Calculate the outflow mass fluxes   !
  !     then correct it to satisfy the    ! 
  !          overall mass balance         !
  !---------------------------------------!
  do m = 1, grid % n_materials
    bulk(m) % mass_out = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW) then
          U % n(c2) = U % n(c1)
          V % n(c2) = V % n(c1)
          W % n(c2) = W % n(c1)
          Flux(s) = DENc(material(c1)) * ( U % n(c2)*grid % sx(s) +  & 
                                           V % n(c2)*grid % sy(s) +  &
                                           W % n(c2)*grid % sz(s) )
          if(material(c1) == m) then
            bulk(m) % mass_out = bulk(m) % mass_out + Flux(s)
          end if
        end if
        if(TypeBC(c2) == CONVECT.and.Flux(s) > 0.0) then
          Flux(s) = DENc(material(c1)) * ( U % n(c2)*grid % sx(s) +  & 
                                           V % n(c2)*grid % sy(s) +  &
                                           W % n(c2)*grid % sz(s) )
          if(material(c1) == m) then
            bulk(m) % mass_out = bulk(m) % mass_out + Flux(s)
          end if
        end if

        Flux(s) = DENc(material(c1)) * ( U % n(c2)*grid % sx(s) +  & 
                                         V % n(c2)*grid % sy(s) +  &
                                         W % n(c2)*grid % sz(s) )
        if(TypeBC(c2) == PRESSURE.and.Flux(s) > 0.0) then
          if(material(c1) == m) then
            bulk(m) % mass_out = bulk(m) % mass_out + Flux(s)
          end if
        endif

      endif
    end do
    call glosum(bulk(m) % mass_out)  ! not checked
  end do

  do m = 1, grid % n_materials
    fac(m) = bulk(m) % mass_in/(bulk(m) % mass_out+TINY)
  end do

  do m = 1, grid % n_materials
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW  .or.  &
           TypeBC(c2) == CONVECT  .or.  &
           TypeBC(c2) == PRESSURE) then
          if(material(c1) == m) then
            U % n(c2) = U % n(c2) * fac(m)
            V % n(c2) = V % n(c2) * fac(m)
            W % n(c2) = W % n(c2) * fac(m)
            Flux(s) = Flux(s)*fac(m) 
            bulk(m) % mass_out = bulk(m) % mass_out + Flux(s)
          endif
        endif
      endif
    end do
  end do

  end subroutine
