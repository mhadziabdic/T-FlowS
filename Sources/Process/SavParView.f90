!======================================================================!
  subroutine SavParView(sub, NCsub, namAut, n0, n1)
!----------------------------------------------------------------------!
! Reads: NAME.gmv and generates NAME.vti Paraview XML output file      !
! ~~~~~~~                                                              ! 
!------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod
  use allp_mod
  use les_mod
  use pro_mod
  use rans_mod
  use Tokenizer_Mod
  use Constants_Pro_Mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer ::  NCsub, sub, n0, n1
  character :: storename*80, namTem*80, namXML*80
  character :: namAut*(*)
!-------------------------------[Locals]-------------------------------!
  integer   :: c, n
  character :: name_out*80, stringadummy*100, nameIn*80
  real,allocatable :: x(:), y(:), z(:)    ! self evident
  integer,allocatable :: connessione(:,:) ! connection
  integer :: celleconnessione
  integer :: NNsub, NCsub_new
  integer :: off_set_connection
  integer :: i
!======================================================================!

  namTem = problem_name
  storename = namAut
!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     reads GMV file      !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  call Name_File(sub, name_out, '.gmv')
  open(9, file=name_out)
  if (this_proc <2) then
    print *, 'Now reading the file: ', name_out
  end if
!---------------!
!     start     !
!---------------!
  call Tokenizer_Mod_Read_Line(9)  !read 'gmvinput ascii' line
!---------------!
!     nodes     !
!---------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1),*) stringadummy
  read(line % tokens(2),*) NNsub
  allocate(x(NNsub)); x=0.
  allocate(y(NNsub)); y=0.
  allocate(z(NNsub)); z=0.
  
  do n=1,NNsub
    call Tokenizer_Mod_Read_Line(9)                           
    read(line % tokens(1),*) x(n)
  end do
  do n=1,NNsub
    call Tokenizer_Mod_Read_Line(9)                           
    read(line % tokens(1),*) y(n)
  end do
  do n=1,NNsub
    call Tokenizer_Mod_Read_Line(9)                           
    read(line % tokens(1),*) z(n)
  end do

!----------------------!
!     cell section     !
!----------------------!

  celleconnessione = 0
  off_set_connection = 0

  call Tokenizer_Mod_Read_Line(9)                           
  read(line % tokens(1),*) stringadummy
  read(line % tokens(2),*) NCsub_new

  if (NCsub_new.ne.NCsub) then 
     print *, 'number of cells read and processed is different, exiting!'
     stop
  end if
  
  allocate(connessione(NCsub,9)); connessione=0

  do n=1,NCsub
    call Tokenizer_Mod_Read_Line(9)                           
    read(line % tokens(1),*) stringadummy
    read(line % tokens(2),*) off_set_connection
    if (n==1) then 
       connessione(n,1) = off_set_connection
    else
       connessione(n,1) = connessione (n-1,1)+ off_set_connection
    end if
    
    celleconnessione = celleconnessione + off_set_connection

    call Tokenizer_Mod_Read_Line(9)
    do c=1,off_set_connection
      read(line % tokens(c),*) connessione(n,c+1)
    end do
  end do

  close(9)

  problem_name = namAut
  call Name_File(sub, namXML, '.vtu')

  open(9, file=namXML)
  if (this_proc <2) then
  write(6, *) 'Now writing the file:', namXML
  end if


  write(9,'(A21)') '<?xml version="1.0"?>'
  write(9,'(A73)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(9,'(A20)') '  <UnstructuredGrid>'

  write(9,*) '    <Piece NumberOfPoints="', NNsub,'" NumberOfCells="', NCsub, '" >' 


  write(9,*) '       <CellData Scalars="scalars" vectors="velocity">'

  !--scalar: pressure
  write(9,*) '        <DataArray type="Float32" Name="pressure" format="ascii">'
  do c=1,NCsub
    write(9,*) P % n(c)
  end do  
  write(9,*) '        </DataArray>'
  
!--vector: velocity
  write(9,'(A160)') '<DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
  do c=1,NCsub
    write(9,*) U % n(c), V % n(c), W % n(c) 
  end do  
  write(9,'(A20)') '        </DataArray>'

  write(9,'(A65)') '<DataArray type="Float32" Name="Y+" format="ascii">'
  do c=1,NCsub
    if(IsNearWall(c)) then 
      write(9,*) sqrt(WallDs(c)*sqrt(U%n(c)**2+V%n(c)**2+W%n(c)**2)/VISc)
    else
      write(9,*) 0.0
    end if
  end do  
  write(9,'(A20)') '        </DataArray>'

  if(HOT == YES) then
    write(9,'(A65)') '        <DataArray type="Float32" Name="T   " format="ascii">'
    do c=1,NCsub
      write(9,*) T % n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
  end if 

  if(SIMULA == LES.or.SIMULA==K_EPS.or.SIMULA==ZETA) then
    write(9,'(A99)') '        <DataArray type="Float32" Name="viscosity ratio" format="ascii">'
    do c=1,NCsub
      write(9,*) VISt(c)/VISc
    end do  
    write(9,'(A20)') '        </DataArray>'
  end if

  if(SIMULA == LES.and.n0 < n1) then
    write(9,'(A65)') '        <DataArray type="Float32" Name="Pmean" format="ascii">'
    do c=1,NCsub
      write(9,*) P % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    
    write(9,'(A99)') '<DataArray type="Float32" Name="mean velocity" NumberOfComponents="3" format="ascii">'
    do c=1,NCsub
      write(9,*) U % mean(c), V % mean(c), W % mean(c) 
    end do  
    write(9,'(A20)') '        </DataArray>'


    write(9,'(A100)') '        <DataArray type="Float32" Name="Turb. Kin. Energy" format="ascii">'
    do c=1,NCsub
      write(9,*) 0.5*(uu % mean(c) - U % mean(c) * U % mean(c) + &
                      vv % mean(c) - V % mean(c) * V % mean(c) + &
                      ww % mean(c) - W % mean(c) * W % mean(c))
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="uu" format="ascii">'
    do c=1,NCsub
      write(9,*) uu % mean(c) - U % mean(c) * U % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="vv" format="ascii">'
      do c=1,NCsub
        write(9,*) vv % mean(c) - V % mean(c) * V % mean(c)
      end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="ww" format="ascii">'
    do c=1,NCsub
      write(9,*) ww % mean(c) - W % mean(c) * W % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="uv" format="ascii">'
    do c=1,NCsub
      write(9,*) uv % mean(c) - U % mean(c) * V % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="uw" format="ascii">'
    do c=1,NCsub
      write(9,*) uw % mean(c) - U % mean(c) * W % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A65)') '        <DataArray type="Float32" Name="vw" format="ascii">'
    do c=1,NCsub
      write(9,*) vw % mean(c) - V % mean(c) * W % mean(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    if(HOT == YES) then
      write(9,'(A65)') '        <DataArray type="Float32" Name="Tmean" format="ascii">'
      do c=1,NCsub
        write(9,*) T % mean(c)
      end do  
      write(9,'(A20)') '        </DataArray>'

      write(9,'(A65)') '        <DataArray type="Float32" Name="TT" format="ascii">'
      do c=1,NCsub
        write(9,*) TT % mean(c) - T % mean(c) * T % mean(c)
      end do  
      write(9,'(A20)') '        </DataArray>'

      write(9,'(A65)') '        <DataArray type="Float32" Name="uT" format="ascii">'
      do c=1,NCsub
        write(9,*) uT % mean(c) - U % mean(c) * T % mean(c)
      end do  
      write(9,'(A20)') '        </DataArray>'

      write(9,'(A65)') '        <DataArray type="Float32" Name="vT" format="ascii">'
      do c=1,NCsub
        write(9,*) vT % mean(c) - V % mean(c) * T % mean(c)
      end do  
      write(9,'(A20)') '        </DataArray>'

      write(9,'(A65)') '        <DataArray type="Float32" Name="wT" format="ascii">'
      do c=1,NCsub
        write(9,*) wT % mean(c) - W % mean(c) * T % mean(c)
      end do  
      write(9,'(A20)') '        </DataArray>'
    end if
  end if

  if(SIMULA == EBM.or.SIMULA==HJ) then
    write(9,'(A99)') '        <DataArray type="Float32" Name="TKE" format="ascii">'
    do c=1,NCsub
      write(9,*) Kin%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A99)') '        <DataArray type="Float32" Name="Eps" format="ascii">'
    do c=1,NCsub
      write(9,*) Eps%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A99)') '        <DataArray type="Float32" Name="uu" format="ascii">'
    do c=1,NCsub
      write(9,*) uu%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A99)') '        <DataArray type="Float32" Name="vv" format="ascii">'
    do c=1,NCsub
      write(9,*) vv%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A99)') '        <DataArray type="Float32" Name="ww" format="ascii">'
    do c=1,NCsub
      write(9,*) ww%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A99)') '        <DataArray type="Float32" Name="uv" format="ascii">'
    do c=1,NCsub
      write(9,*) uv%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A99)') '        <DataArray type="Float32" Name="uw" format="ascii">'
    do c=1,NCsub
      write(9,*) uw%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
    write(9,'(A99)') '        <DataArray type="Float32" Name="vw" format="ascii">'
    do c=1,NCsub
      write(9,*) vw%n(c)
    end do  
    write(9,'(A20)') '        </DataArray>'
  end if    
  if(SIMULA == K_EPS.or.SIMULA==ZETA) then
    write(9,'(A99)') '        <DataArray type="Float32" Name="TKE" format="ascii">'
    do c=1,NCsub
      write(9,*) Kin%n(c) 
    end do  
    write(9,'(A20)') '        </DataArray>'

    write(9,'(A99)') '        <DataArray type="Float32" Name="EPS" format="ascii">'
    do c=1,NCsub
      write(9,*) Eps%n(c) 
    end do  
    write(9,'(A20)') '        </DataArray>'

    if(BUOY == YES) then
      write(9,'(A99)') '        <DataArray type="Float32" Name="tt" format="ascii">'
      do c=1,NCsub
        write(9,*) tt%n(c) 
      end do  
      write(9,'(A20)') '        </DataArray>'
      write(9,'(A99)') '        <DataArray type="Float32" Name="ut" format="ascii">'
      do c=1,NCsub
        write(9,*) ut%n(c) 
      end do  
      write(9,'(A20)') '        </DataArray>'
      write(9,'(A99)') '        <DataArray type="Float32" Name="vt" format="ascii">'
      do c=1,NCsub
        write(9,*) vt%n(c) 
      end do  
      write(9,'(A20)') '        </DataArray>'
      write(9,'(A99)') '        <DataArray type="Float32" Name="wt" format="ascii">'
      do c=1,NCsub
        write(9,*) wt%n(c) 
      end do  
      write(9,'(A20)') '        </DataArray>'
    end if
  end if

  call wait

  write(9,'(A17)') '      </CellData>'
  
  write(9,'(A14)') '      <Points>'
  write(9,*) '         <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
                             do c=1,NNsub
                                 write(9,*) x(c), y(c), z(c)
                             end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A15)') '      </Points>'
  write(9,'(A13)') '      <Cells>'
  write(9,'(A71)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do n=1,NCsub
                                   !hexa
                                   if (connessione(n,9) /= 0) then
                                     do i=2,9
                                            write(9,*) connessione(n,i)-1
                                     end do
                                   else
                                       !prism
                                       if ((connessione(n,8) == 0).and.(connessione(n,7) /= 0)) then
                                         write(9,*) connessione(n,2)-1, connessione(n,4)-1, & 
                                         connessione(n,3)-1, connessione(n,5)-1, connessione(n,7)-1, connessione(n,6)-1
                                       else
                                          if ((connessione(n,7) == 0).and.(connessione(n,6) /= 0)) then
                                              write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                                              connessione(n,3)-1, connessione(n,6)-1, connessione(n,2)-1
                                          else
                                              write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                                              connessione(n,3)-1, connessione(n,2)-1
                                          end if
                                       end if
                                   end if

                              end do 
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A62)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
                              do n=1,NCsub
                                 write(9,*) connessione(n,1)
                              end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A60)') '        <DataArray type="UInt8" Name="types" format="ascii">'
                              do n=1,NCsub
                                 if ((connessione(n,6)==0)) then !thetra
                                      write(9,*) '10'
                                 else 
                                 if ((connessione(n,7)==0)) then
                                      write(9,*) '14'          !pyr
                                        else
                                        if ((connessione(n,8)==0)) then
                                            write(9,*) '13'      !prism
                                        else
                                            write(9,*) '12'      !hexa
                                        end if
                                     end if
                                 end if
                              end do
  write(9,'(A20)') '        </DataArray>'
  write(9,'(A14)') '      </Cells>'
  write(9,'(A14)') '</Piece>'
  write(9,'(A19)') '</UnstructuredGrid>'
  write(9,'(A14)') '</VTKFile>'

! write down the master file for parallel jobs
! actually it writes it down also for sequential process but it's not
! necessary for paraview

  if (n_proc > 1) then
    if(this_proc < 2) then
    problem_name = storename
    call Name_File(0, problem_name, '.pvtu')
    open(112, file=problem_name)
    write(112,'(A21)') '<?xml version="1.0"?>'
    write(112,'(A74)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(112,*) '<PUnstructuredGrid GhostLevel="0">'
    write(112,*) '       <PCellData Scalars="scalars" vectors="velocity">'
    write(112,*) '        <PDataArray type="Float32" Name="pressure"/>'
    write(112,*) '        <PDataArray type="Float32" Name="velocity" NumberOfComponents="3"/>'
    write(112,*) '        <PDataArray type="Float32" Name="Y+"/>'
    if(HOT == YES) write(112,*) '        <PDataArray type="Float32" Name="T"/>'
    if(SIMULA == LES) then
      write(112,*) '        <PDataArray type="Float32" Name="Viscosity ratio"/>'
      if(n0<n1) then
        write(112,*) '        <PDataArray type="Float32" Name="mean pressure"/>'
        write(112,'(A80)') '<PDataArray type="Float32" Name="mean velocity" NumberOfComponents="3"/>'
        write(112,*) '        <PDataArray type="Float32" Name="Turb. Kin. Energy"/>'
        write(112,*) '        <PDataArray type="Float32" Name="uu"/>'
        write(112,*) '        <PDataArray type="Float32" Name="vv"/>'
        write(112,*) '        <PDataArray type="Float32" Name="ww"/>'
        write(112,*) '        <PDataArray type="Float32" Name="uv"/>'
        write(112,*) '        <PDataArray type="Float32" Name="uw"/>'
        write(112,*) '        <PDataArray type="Float32" Name="vw"/>'
        if(HOT == YES) then
          write(112,'(A65)') '        <PDataArray type="Float32" Name="Tmean"/>'
          write(112,'(A65)') '        <PDataArray type="Float32" Name="TT"/>'
          write(112,'(A65)') '        <PDataArray type="Float32" Name="Tu"/>'
          write(112,'(A65)') '        <PDataArray type="Float32" Name="Tv"/>'
          write(112,'(A65)') '        <PDataArray type="Float32" Name="Tw"/>'
        end if
      end if
    end if
    if(SIMULA == K_EPS.or.SIMULA==ZETA) then
      write(112,*) '        <PDataArray type="Float32" Name="TKE"/>'
      write(112,*) '        <PDataArray type="Float32" Name="EPS"/>'
      write(112,*) '        <PDataArray type="Float32" Name="uw"/>'
    end if 
    if(SIMULA == EBM.or.SIMULA==HJ) then
      write(112,*) '        <PDataArray type="Float32" Name="TKE"/>'
      write(112,*) '        <PDataArray type="Float32" Name="EPS"/>'
      write(112,*) '        <PDataArray type="Float32" Name="uu"/>'
      write(112,*) '        <PDataArray type="Float32" Name="vv"/>'
      write(112,*) '        <PDataArray type="Float32" Name="ww"/>'
      write(112,*) '        <PDataArray type="Float32" Name="uv"/>'
      write(112,*) '        <PDataArray type="Float32" Name="uw"/>'
      write(112,*) '        <PDataArray type="Float32" Name="vw"/>'
    end if
  end if
 
  write(112,*) '       </PCellData>'
  write(112,*) '       <PPoints>'
  write(112,*) '         <PDataArray type="Float32" NumberOfComponents="3"/>'
  write(112,*) '       </PPoints>'
  problem_name = namAut
  do i=1,n_proc
    call Name_File(i, nameIn, '.vtu"/>')
    write(112,'(A)',advance="no") '<Piece Source="'
    write(112,'(A)') nameIn
  end do
     
  write(112,*) '</PUnstructuredGrid>'
  write(112,*) '</VTKFile>'
  close(112)
  end if !(this_proc <2)

  close(9)
  problem_name = namTem

  call wait

  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(connessione)

  return

  end subroutine
