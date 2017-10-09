!======================================================================!
  subroutine DivLoa
!----------------------------------------------------------------------!
! Reads:  NAME.cns, NAME.gmv                                           !
! ~~~~~~                                                               !                                                                           
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer           :: c, n, s
  character(len=80) :: dum_s, name_in
!======================================================================!

  write(*,'(A41)') '# Input problem name: (without extension)'
  call ReadC(5,inp,tn,ts,te)  
  read(inp, '(A80)')  name

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!       Read the file with the      !
!     connections between cells     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.cns'

  open(9, FILE=name_in,FORM='UNFORMATTED')
  write(*,'(A24,A)') '# Now reading the file: ', name_in

!///// number of cells
  read(9) NC
  read(9) NbC
  read(9) NS
  read(9) NSsh
  read(9) Nmat

  allocate (material(-NbC:NC))
  read(9) (material(c), c=1,NC)
  read(9) (material(c), c=-1,-NBC,-1)

!///// sides
  allocate (SideC(0:2,NS+NSsh))
  read(9) (SideC(0,s), s=1,NS)
  read(9) (SideC(1,s), s=1,NS)
  read(9) (SideC(2,s), s=1,NS)

!///// boundary cells
  allocate (bcmark(-NbC-1:-1)); bcmark=0
  allocate (CopyC(-NbC:-1));    CopyC=0
  read(9) (bcmark(c), c=-1,-NbC, -1)
  read(9) (CopyC(c), c=-1,-NbC, -1)
 
  read(9) Ncopy
  write(*,*) Ncopy
  allocate (CopyS(2,Ncopy))
  read(9) (CopyS(1,s), s=1,Ncopy)
  read(9) (CopyS(2,s), s=1,Ncopy)

  close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Read the GMV file with node      !
!     coordinates and volume nodes     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.gmv'
  write(*,*) '# Now reading the file:', name_in
  open(9, FILE=name_in)

!>>>>> read the number of nodes
  call ReadC(9,inp,tn,ts,te)
  call ReadC(9,inp,tn,ts,te)
  read(inp(ts(2):te(2)),*)   NN   ! number of nodes

!////////////////////////!
!     Alocate memory     ! 
!////////////////////////!

  write(6,'(A25)')       '# Allocating memory for: ' 
  write(6,'(A1,I8,A6)')  '#', NN,  ' nodes' 
  write(6,'(A1,I8,A6)')  '#', NC,  ' cells' 
  write(6,'(A1,I8,A15)') '#', NbC, ' boundary cells'         
  write(6,'(A1,I8,A11)') '#', NS,  ' cell faces' 

!---- variables defined in all_mod.h90:
  allocate (xc(-NbC:NC)); xc=0.0
  allocate (yc(-NbC:NC)); yc=0.0
  allocate (zc(-NbC:NC)); zc=0.0
  allocate (Sx(NS)); Sx=0.0
  allocate (Sy(NS)); Sy=0.0
  allocate (Sz(NS)); Sz=0.0
  allocate (Dx(NS+NSsh)); Dx=0.0
  allocate (Dy(NS+NSsh)); Dy=0.0
  allocate (Dz(NS+NSsh)); Dz=0.0
  allocate (xsp(NS)); xsp=0.0
  allocate (ysp(NS)); ysp=0.0
  allocate (zsp(NS)); zsp=0.0
  allocate (volume(-NbC:NC)); volume=0.0
  allocate (delta(-NbC:NC));  delta=0.0
  allocate (WallDs(NS)); WallDs=0.0
  allocate (f(NS)); f=0.0

!---- variables declared in gen_mod.h90:
  allocate (NewC(-NbC-1:NC)); NewC=0 
  allocate (NewS(NS));        NewS=0

  allocate (x_node(NN)); x_node=0
  allocate (y_node(NN)); y_node=0
  allocate (z_node(NN)); z_node=0

  allocate (NewN(NN)); NewN=0 

  allocate (CellN(NC,0:8)); CellN=0
  allocate (SideN(NS+NSsh,0:4)); SideN=0

!---- variables declared in div.h90:
  allocate (ix(-NbC:NC));  ix=0
  allocate (iy(-NbC:NC));  iy=0
  allocate (iz(-NbC:NC));  iz=0
  allocate (iin(-NbC:NC)); iin=0
  allocate (criter(NC));   criter=0

!---- variables declared in par_mod.h90:
  allocate (proces(NC)); proces=0
  allocate (BuSeIn(NS)); BuSeIn=0
  allocate (BuReIn(NS)); BuReIn=0
  allocate (BufPos(NS)); BufPos=0

  write(6,'(A26)') '# Allocation successfull !'

!>>>>> read node coordinates
  do n=1,NN
    read(9,*) x_node(n)
  end do
  do n=1,NN
    read(9,*) y_node(n)
  end do
  do n=1,NN
    read(9,*) z_node(n)
  end do

!>>>>> read cell nodes 
  call ReadC(9,inp,tn,ts,te)  ! cells, number of cells
  do c=1,NC  !->>> ovo bi se moglo napravit neformatirano
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) dum_s
    if(dum_s  ==  'hex') then
      CellN(c,0) = 8
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) &
           CellN(c,1), CellN(c,2), CellN(c,4), CellN(c,3),          &
           CellN(c,5), CellN(c,6), CellN(c,8), CellN(c,7) 
    else if(dum_s  ==  'prism') then
      CellN(c,0) = 6
      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                             &
         CellN(c,1), CellN(c,2), CellN(c,3),  &
         CellN(c,4), CellN(c,5), CellN(c,6)
    else if(dum_s  ==  'tet') then
      CellN(c,0) = 4
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) &
         CellN(c,1), CellN(c,2), CellN(c,3), CellN(c,4)
    else if(dum_s  ==  'pyramid') then
      CellN(c,0) = 5
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) &
         CellN(c,5), CellN(c,1), CellN(c,2), CellN(c,4), CellN(c,3)
    else
      write(*,*) 'Unsupported cell type: ', dum_s
      write(*,*) 'Exiting'
      stop
    end if
  end do

!>>>>> read cell materials
  call ReadC(9,inp,tn,ts,te)  ! materials, number of materials 
  read(inp(ts(2):te(2)),*) Nmat
  do n=1,Nmat
    call ReadC(9,inp,tn,ts,te)  
  end do 
  do c=1,NC
    read(9,*) material(c)
  end do 

  close(9)

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something ! 

  end subroutine DivLoa