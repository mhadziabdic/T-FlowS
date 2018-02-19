!======================================================================!
  subroutine CalcMn_Cylind(grid, n0, n1)   
!----------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  type(Grid_Type) :: grid
  integer         :: n0, n1
!-------------------------------[Locals]-------------------------------!
  integer           :: c, n
  real              :: Urad_mean, Utan_mean, R, Urad, Utan
!======================================================================!

  call Control_Mod_Heat_Transfer(verbose = .true.)

  n=n1-n0

  if(n  > -1) then
    do c=-grid % n_bnd_cells,grid % n_cells
!-----------------------!
!      mean values      !
!-----------------------!
      R           = (grid % xc(c)*grid % xc(c) + grid % yc(c)*grid % yc(c))**0.5 + tiny
      U % mean(c) = ( U % mean(c) * (1.*n) + U % n(c) ) / (1.*(n+1))
      V % mean(c) = ( V % mean(c) * (1.*n) + V % n(c) ) / (1.*(n+1))
      W % mean(c) = ( W % mean(c) * (1.*n) + W % n(c) ) / (1.*(n+1))
      P % mean(c) = ( P % mean(c) * (1.*n) + P % n(c) ) / (1.*(n+1))

      Urad        = (U % n(c) * grid % xc(c) / R  + V % n(c) * grid % yc(c) / R)
      Utan        = (-U % n(c) * grid % yc(c) / R  + V % n(c) * grid % xc(c) / R)

!------------------------------!
!      fluctuating values      !
!------------------------------!
      uu % mean(c) = ( uu % mean(c)*(1.*n) + Urad * Urad ) &
                   / (1.*(n+1))
      vv % mean(c) = ( vv % mean(c)*(1.*n) + Utan * Utan ) &
                   / (1.*(n+1))
      ww % mean(c) = ( ww % mean(c)*(1.*n) + W % n(c) * W % n(c) ) &
                   / (1.*(n+1))
      uv % mean(c) = ( uv % mean(c)*(1.*n) + Urad * Utan ) &
                   / (1.*(n+1))
      uw % mean(c) = ( uw % mean(c)*(1.*n) + Urad * W % n(c) ) &
                   / (1.*(n+1))
      vw % mean(c) = ( vw % mean(c)*(1.*n) + Utan * W % n(c) ) &
                   / (1.*(n+1))


      vis_t_mean(c) = ( vis_t_mean(c)*(1.*n) + vis_t(c) ) & 
                   / (1.*(n+1))

      if(heat_transfer == YES) then
        T % mean(c) = ( T % mean(c) * (1.*n) + T % n(c) ) / (1.*(n+1))

        TT % mean(c) = ( TT % mean(c)*(1.*n) + T % n(c) * T % n(c) ) & 
                     / (1.*(n+1))
        uT % mean(c) = ( uT % mean(c)*(1.*n) + Urad * T % n(c) ) &
                     / (1.*(n+1))
        vT % mean(c) = ( vT % mean(c)*(1.*n) + Utan * T % n(c) ) &
                     / (1.*(n+1))
        wT % mean(c) = ( wT % mean(c)*(1.*n) + w % n(c) * T % n(c) ) &
                     / (1.*(n+1))
      end if
    end do 
  end if

  if(n  > -1) then
    do c=-grid % n_bnd_cells,grid % n_cells
      R           = (grid % xc(c)*grid % xc(c) + grid % yc(c)*grid % yc(c))**0.5 + tiny
      Urad_mean   = (U % mean(c) * grid % xc(c) / R  + V % mean(c) * grid % yc(c) / R)
      Utan_mean   = (-U % mean(c) * grid % yc(c) / R  + V % mean(c) * grid % xc(c) / R)
      Urad        = (U % n(c) * grid % xc(c) / R  + V % n(c) * grid % yc(c) / R)
      Utan        = (-U % n(c) * grid % yc(c) / R  + V % n(c) * grid % xc(c) / R)

      uuu % mean(c) = ( uuu % mean(c) * (1.*n) +   &
      (Urad - Urad_mean)**3.0 ) / (1.*(n+1))

      uuv % mean(c) = ( uuv % mean(c) * (1.*n) +   &
      (Urad - Urad_mean)**2.0*(Utan - Utan_mean)  )  / (1.*(n+1))

      uuw % mean(c) = ( uuw % mean(c) * (1.*n) +   &
      (Urad - Urad_mean)**2.0*(W % n(c) - W % mean(c) )  )  / (1.*(n+1))

      vvu % mean(c) = ( vvu % mean(c) * (1.*n) +   &
      (Utan - Utan_mean)**2.0*(Urad - Urad_mean)  )  / (1.*(n+1))

      vvv % mean(c) = ( vvv % mean(c) * (1.*n) +   &
      (Utan - Utan_mean)**3.0 )  / (1.*(n+1))

      vvw % mean(c) = ( vvw % mean(c) * (1.*n) +   &
      (Utan - Utan_mean)**2.0*(W % n(c) - W % mean(c))  )  / (1.*(n+1))

      wwu % mean(c) = ( wwu % mean(c) * (1.*n) +   &
      (W % n(c) - W % mean(c))**2.0*(Urad - Urad_mean) )  / (1.*(n+1))

      wwv % mean(c) = ( wwv % mean(c) * (1.*n) +   &
      (W % n(c) - W % mean(c))**2.0*(Utan - Utan_mean)  )  / (1.*(n+1))

      www % mean(c) = ( www % mean(c) * (1.*n) +   &
      (W % n(c) - W % mean(c))**3.0  )  / (1.*(n+1))

      uwu % mean(c) = ( uwu % mean(c) * (1.*n) +   &
      (Urad - Urad_mean)**2.0*(W % n(c) - W % mean(c))  )  / (1.*(n+1))

      uwv % mean(c) = ( uwv % mean(c) * (1.*n) +   &
      (Urad - Urad_mean)*(W % n(c) - W % mean(c))*(Utan - Utan_mean))  / (1.*(n+1))

      uww % mean(c) = ( uww % mean(c) * (1.*n) +   &
      (W % n(c) - W % mean(c))**2.0*(Urad - Urad_mean)  )  / (1.*(n+1))

    end do
  end if

  end subroutine
