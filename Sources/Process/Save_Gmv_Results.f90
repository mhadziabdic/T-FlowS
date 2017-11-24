!==============================================================================!
  subroutine Save_Gmv_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   Writes: NAME.r.gmv                                                         !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod, only: this_proc
  use rans_mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c
  character(len=80)   :: name_out, answer, store_name
  character, optional :: name_save*(*)
!==============================================================================!

  ! Store the name
  store_name = name     

  if(PRESENT(name_save)) then
    write(*,*) name_save
    name = name_save  
  else
    if(this_proc  < 2)  &
      write(*,*) '# Input result file name [skip cancels]:'
    call Tokenizer_Mod_Read_Line(CMN_FILE)  
    read(line % tokens(1),'(A80)')  name
    answer=name
    call To_Upper_Case(answer) 
    if(answer == 'SKIP') then
      name = store_name  
      return
    end if 
  end if

  call wait 

  !------------------------------!
  !                              !
  !   Create .gmv results file   !
  !                              !
  !------------------------------!
  call Name_File(this_proc, name_out, '.r.gmv', len_trim('.r.gmv'))
  open(9, file=name_out)
  write(*,*) '# Now creating the file:', trim(name_out)

  !--------------!
  !   Velocity   !
  !--------------!
  write(9, *) 'velocity 0'  !  0 is for cell data
  write(9,'(1F16.6)') (U % n(c),c=1,grid % n_cells)
  write(9,'(1F16.6)') (V % n(c),c=1,grid % n_cells)
  write(9,'(1F16.6)') (W % n(c),c=1,grid % n_cells)

  !---------------------!
  !   Other variables   !
  !---------------------!
  write(9, *) 'variable' 

  !--------------!
  !   Pressure   !
  !--------------!
  write(9, *) 'P 0'
  write(9,'(1F16.6)') (P % n(c),c=1,grid % n_cells)

  !--------------------------!
  !   Pressure corrections   !
  !--------------------------!
  write(9, *) 'PP 0'
  write(9,'(1F16.6)') (PP % n(c),c=1,grid % n_cells)

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(HOT==YES) then
    write(9, *) 'T 0'
    write(9,'(1F16.6)') (T % n(c),c=1,grid % n_cells)
  end if

  !-----------!
  !   K_EPS   !
  !-----------!
  if(SIMULA == K_EPS.or.SIMULA == HYB_PITM) then 
    write(9, *) 'k 0'
    write(9,'(1F16.6)') (Kin % n(c),c=1,grid % n_cells)
    write(9, *) 'eps 0'
    write(9,'(1F16.6)') (Eps % n(c),c=1,grid % n_cells)
    write(9, *) 'VISt 0'
    write(9,'(1F16.6)') (VISt(c),c=1,grid % n_cells)
    write(9, *) 'Pk 0'
    write(9,'(1F16.6)') (Pk(c),c=1,grid % n_cells)
    write(9, *) 'Yplus 0'
    write(9,'(1F16.6)') (Ynd(c),c=1,grid % n_cells)
    write(9, *) 'Uf 0'
    write(9,'(1F16.6)') (Uf(c),c=1,grid % n_cells)
  end if 

  !--------------!
  !   K_EPS_VV   !
  !--------------!
  if(SIMULA == K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then 
    write(9, *) 'k 0'
    write(9,'(1F16.6)') (Kin % n(c),c=1,grid % n_cells)
    write(9, *) 'eps 0'
    write(9,'(1F16.6)') (Eps % n(c),c=1,grid % n_cells)
    write(9, *) 'vv 0'
    write(9,'(1F16.6)') (v_2 % n(c),c=1,grid % n_cells)
    write(9, *) 'f22 0'
    write(9,'(1F16.6)') (f22 % n(c),c=1,grid % n_cells)
    write(9, *) 'Pk 0'
    write(9,'(1F16.6)') ( Pk(c),c=1,grid % n_cells)
    write(9, *) 'VISt 0'
    write(9,'(1F16.6)') (VISt(c),c=1,grid % n_cells)
  end if 

  !-------------!
  !   SPA_ALL   !
  !-------------!
  if(SIMULA == SPA_ALL) then
    write(9, *) 'VIS 0'
    write(9,'(1F16.6)') (VIS % n(c),c=1,grid % n_cells)
    write(9, *) 'Vort 0'
    write(9,'(1F16.6)') (Vort(c),c=1,grid % n_cells)
    write(9, *) 'VISt 0'
    write(9,'(1F16.6)') (VISt(c),c=1,grid % n_cells)
    write(9, *) 'WallDs 0'
    write(9,'(1F16.6)') (WallDs(c),c=1,grid % n_cells)
    write(9, *) 'delta 0'
    write(9,'(1F16.6)') (delta(c),c=1,grid % n_cells)
  end if

  !-------------!
  !   DES_SPA   !  
  !-------------!
  if(SIMULA == DES_SPA) then

    ! Mean velocities 
    write(9, *) 'Umean 0'
    write(9,'(1F16.6)') (U % mean(c),c=1,grid % n_cells)
    write(9, *) 'Vmean 0'
    write(9,'(1F16.6)') (V % mean(c),c=1,grid % n_cells)
    write(9, *) 'Wmean 0'
    write(9,'(1F16.6)') (W % mean(c),c=1,grid % n_cells)
    if(HOT == YES) then
      write(9, *) 'Tmean 0'
      write(9,'(1F16.6)') (T % mean(c),c=1,grid % n_cells)
    end if

    ! Velocity fluctuations
    write(9, *) 'uu 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UU % mean(c)-U % mean(c)*U % mean(c),c=1,grid % n_cells)
    write(9, *) 'vv 0'  !  0 is for cell data
    write(9,'(1F16.6)')(VV % mean(c)-V % mean(c)*V % mean(c),c=1,grid % n_cells)
    write(9, *) 'ww 0'  !  0 is for cell data
    write(9,'(1F16.6)')(WW % mean(c)-W % mean(c)*W % mean(c),c=1,grid % n_cells)
    write(9, *) 'uv 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UV % mean(c)-U % mean(c)*V % mean(c),c=1,grid % n_cells)
    write(9, *) 'uw 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UW % mean(c)-U % mean(c)*W % mean(c),c=1,grid % n_cells)
    write(9, *) 'vw 0'  !  0 is for cell data
    write(9,'(1F16.6)')(VW % mean(c)-V % mean(c)*W % mean(c),c=1,grid % n_cells)

    write(9, *) 'VIS 0'
    write(9,'(1F16.6)') (VIS % n(c),c=1,grid % n_cells)
    write(9, *) 'Vort 0'
    write(9,'(1F16.6)') (Vort(c),c=1,grid % n_cells)
    write(9, *) 'VISt 0'
    write(9,'(1F16.6)') (VISt(c),c=1,grid % n_cells)
    write(9, *) 'WallDs 0'
    write(9,'(1F16.6)') (WallDs(c),c=1,grid % n_cells)
    write(9, *) 'delta 0'
    write(9,'(1F16.6)') (delta(c),c=1,grid % n_cells)
  end if

  !---------!
  !   LES   !
  !---------!
  if(SIMULA == LES) then 

    ! Mean velocities 
    write(9, *) 'Umean 0'
    write(9,'(1F16.6)') (U % mean(c),c=1,grid % n_cells)
    write(9, *) 'Vmean 0'
    write(9,'(1F16.6)') (V % mean(c),c=1,grid % n_cells)
    write(9, *) 'Wmean 0'
    write(9,'(1F16.6)') (W % mean(c),c=1,grid % n_cells)
    if(HOT == YES) then
      write(9, *) 'Tmean 0'
      write(9,'(1F16.6)') (T % mean(c),c=1,grid % n_cells)
    end if

    ! Velocity fluctuations
    write(9, *) 'uu 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UU % mean(c)-U % mean(c)*U % mean(c),c=1,grid % n_cells)
    write(9, *) 'vv 0'  !  0 is for cell data
    write(9,'(1F16.6)')(VV % mean(c)-V % mean(c)*V % mean(c),c=1,grid % n_cells)
    write(9, *) 'ww 0'  !  0 is for cell data
    write(9,'(1F16.6)')(WW % mean(c)-W % mean(c)*W % mean(c),c=1,grid % n_cells)
    write(9, *) 'uv 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UV % mean(c)-U % mean(c)*V % mean(c),c=1,grid % n_cells)
    write(9, *) 'uw 0'  !  0 is for cell data
    write(9,'(1F16.6)')(UW % mean(c)-U % mean(c)*W % mean(c),c=1,grid % n_cells)
    write(9, *) 'vw 0'  !  0 is for cell data
    write(9,'(1F16.6)')(VW % mean(c)-V % mean(c)*W % mean(c),c=1,grid % n_cells)

    ! Turbulent viscosity
    write(9, *) 'muT 0'
    write(9,'(1F16.6)') (VISt(c),c=1,grid % n_cells)

    ! Turbulent viscosity
    write(9, *) 'ShearMean 0'
    write(9,'(1F16.6)') (ShearMean(c),c=1,grid % n_cells)

    ! Wall distance            
    write(9, *) 'wall 0'
    write(9,'(1F16.6)') (WallDs(c),c=1,grid % n_cells)
  end if

  write(9, *) 'endvars' 

  write(9, *) 'endgmv'

  ! This is important for parallel version
  write(9,'(I8)') grid % n_cells    

  close(9)

  ! Restore the name
  name = store_name

  end subroutine
