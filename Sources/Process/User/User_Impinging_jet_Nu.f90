!=======================================================================!
  subroutine User_Impinging_jet_Nu(grid)     
!-----------------------------------------------------------------------!
! The subroutine creates ascii file with Nusselt number that is averaged!
! in the azimuthal direction.                                           !    
!-----------------------------------------------------------------------!
  use all_mod
  use allp_mod
  use les_mod
  use pro_mod
  use par_mod
  use rans_mod
  use Grid_Mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  type(Grid_Type) :: grid
!------------------------------[Calling]-------------------------------!
  interface
    logical function Approx(A,B,tol)
      real           :: A,B
      real, optional :: tol
    end function Approx
  end interface
!-------------------------------[Locals]-------------------------------!
  integer             :: Nprob, pl, i, count, s, c1, c2
  real,allocatable    :: z_p(:), Ump(:), Vmp(:), Wmp(:), &
                               Tmp(:),                 &
                               ind(:),       &
                               var_1(:), var_2(:), Rad_mp(:),     &
                               var_3(:), Rad_1(:), &
                               var_4(:), var_5(:), var_6(:)
  integer,allocatable :: Np(:), Ncount(:), Ncount2(:)
  real                :: R, Rad_2
  logical :: there
!======================================================================!

  INQUIRE( file='rad_coordinate.dat', EXIST=THERE ) 
  if(.NOT.THERE) then
    if(this_proc < 2) write(*,*) "==================================================================="
    if(this_proc < 2) write(*,*) "In order to extract Nusselt number profile                         "
    if(this_proc < 2) write(*,*) "an ascii file with cell-faces coordinates has to be read.          "
    if(this_proc < 2) write(*,*) "The name of the file is rad_coordinate.dat.                        "
    if(this_proc < 2) write(*,*) "The file format should be as follows:                              "
    if(this_proc < 2) write(*,*) "10  ! number of cells + 1                                          "
    if(this_proc < 2) write(*,*) "0.0                                                                "
    if(this_proc < 2) write(*,*) "0.1                                                                "
    if(this_proc < 2) write(*,*) "0.2                                                                "
    if(this_proc < 2) write(*,*) "...                                                                "
    if(this_proc < 2) write(*,*) "==================================================================="
    return
  end if

  open(9, file='rad_coordinate.dat')
!---- write the number of searching intervals 
  read(9,*) Nprob
  allocate(z_p(Nprob*2))
  allocate(ind(Nprob*2))

!---- read the intervals positions
  do pl=1,Nprob
    read(9,*) ind(pl), z_p(pl) 
  end do
  close(9)

  allocate(Np(Nprob));     Np=0
  allocate(Ump(Nprob));    Ump=0.0
  allocate(Vmp(Nprob));    Vmp=0.0
  allocate(Wmp(Nprob));    Wmp=0.0
  allocate(Rad_mp(Nprob));  Rad_mp=0.0
  allocate(var_1(Nprob));  var_1=0.0
  allocate(var_2(Nprob));  var_2=0.0
  allocate(var_3(Nprob));  var_3=0.0
  allocate(var_4(Nprob));  var_4=0.0
  allocate(var_5(Nprob));  var_5=0.0
  allocate(var_6(Nprob));  var_6=0.0

  allocate(Rad_1(Nprob));  Rad_1=0.0
  allocate(Ncount(Nprob)); Ncount=0
  allocate(Ncount2(Nprob)); Ncount2=0

  count = 0

  if(HOT==YES) then
    allocate(Tmp(Nprob));   Tmp=0.0
  end if  

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!
  do i = 1, Nprob
    Rad_1(i) = abs(z_p(i))
  end do

  do i = 1, Nprob-1
    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)
      if(c2 < 0) then
        if(TypeBC(c2).eq.WALLFL) then
          Rad_2 = (grid % xc(c1)*grid % xc(c1) + &
          grid % yc(c1)*grid % yc(c1))**0.5 + tiny
          if(Rad_2 < Rad_1(i+1).and.Rad_2 > Rad_1(i).and.&
            grid % zc(c1) < 0.5) then
            R        = (grid % xc(c1)*grid % xc(c1) + &
                        grid % yc(c1)*grid % yc(c1))**0.5 + tiny
            Rad_mp(i)= Rad_mp(i) + (grid % xc(c1)*grid % xc(c1) +&
                                    grid % yc(c1)*grid % yc(c1))**0.5
            Ump(i)   = Ump(i) + U % n(c1) * grid % xc(c1) / R  + &
                                V % n(c1) * grid % yc(c1) / R
            Vmp(i)   = Vmp(i) + (-U % n(c1) * grid % yc(c1) / R  + &
                                  V % n(c1) * grid % xc(c1) / R)
            Wmp(i)   = Wmp(i) + W % n(c1)
            Tmp(i)   = Tmp(i) + T % n(c2) 
            var_1(i) = var_1(i) + grid % zc(c1)
            var_2(i) = var_2(i) + TauWall(c1)**0.5 
            var_3(i) = var_3(i) + (Cmu**0.25*Kin%n(c1)**0.5) 
            var_4(i) = var_4(i) + Kin %n(c1)
            var_5(i) = var_5(i) + Eps%n(c1)
            var_6(i) = var_6(i) + T % q(c2)
            Ncount(i)= Ncount(i) + 1
          end if
        end if
      end if
    end do
  end do
!---- average over all processors
  do pl=1, Nprob
    call IGlSum(Ncount(pl))
    call IGlSum(Ncount2(pl))

    call GloSum(Ump(pl))
    call GloSum(Vmp(pl))
    call GloSum(Wmp(pl))

    call GloSum(var_1(pl))
    call GloSum(var_2(pl))
    call GloSum(var_3(pl))
    call GloSum(var_4(pl))
    call GloSum(var_5(pl))
    call GloSum(var_6(pl))

    call GloSum(Rad_mp(pl))
    call GloSum(Tmp(pl))

    count =  count + Ncount(pl)
  end do

  do i = 1, Nprob
    if(Ncount(i) /= 0) then
      Wmp(i)    =  Wmp(i)/Ncount(i)
      Vmp(i)    =  Vmp(i)/Ncount(i)
      Ump(i)    =  Ump(i)/Ncount(i)
      Tmp(i)    =  Tmp(i)/Ncount(i)
      var_1(i)  =  var_1(i)/Ncount(i)
      var_2(i)  =  var_2(i)/Ncount(i)
      var_3(i)  =  var_3(i)/Ncount(i)
      var_4(i)  =  var_4(i)/Ncount(i)
      var_5(i)  =  var_5(i)/Ncount(i)
      var_6(i)  =  var_6(i)/Ncount(i)
      Rad_mp(i) =  Rad_mp(i)/Ncount(i)
    end if
  end do
  call wait

  open(3,file='Nusselt_Utau.dat')
  write(3,*) '# Xrad, Nu, Ufriction, yPlus, Temp, Numb of points '
  do i = 1, Nprob
    if(Ncount(i) /= 0) then
      write(3,'(5E15.7,I6)') Rad_mp(i)/2.0, 2.0*var_6(i)/(CONc(material(1))*(Tmp(i)-20.0)), &
                          var_2(i), var_2(i)*var_1(i)/VISc, Tmp(i),  &
                          Ncount(i)
    end if
  end do
  close(3)

  deallocate(Np)
  deallocate(Ump)
  deallocate(Vmp)
  deallocate(Wmp)

  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  deallocate(var_4)
  deallocate(var_5)
  deallocate(var_6)
  deallocate(Rad_mp)
  deallocate(Rad_1)
  deallocate(Ncount)
  deallocate(Ncount2)

  if(HOT==YES) then
    deallocate(Tmp)
  end if

  if(this_proc < 2) write(*,*)'Finished with User_Impinging_jet_Nu'

  end subroutine User_Impinging_jet_Nu
