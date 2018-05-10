!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for no MPI functionality.                                           !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
! include 'mpif.h'
!==============================================================================!

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 8
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! These names are ugly but mean number of buffer boundaries start and end
  integer, allocatable :: nbb_s(:), nbb_e(:)
  integer, allocatable :: buffer_index(:)

  ! Variables which follow are for backup saving to single file
  integer                      :: nc_s   ! number of cells in subdomain
  integer                      :: nb_s   ! number of bundary cells in subdoima
  integer                      :: nf_s   ! number of faces in subdomain
  integer                      :: nbf_s  ! number of buffer faces in subdomain
  integer                      :: nc_t   ! total number of cells 
  integer                      :: nb_t   ! total number of bundary cells
  integer                      :: nf_t   ! total number of faces
  integer(kind=4), allocatable :: cell_map(:)
  integer(kind=4), allocatable :: bnd_cell_map(:)
  integer(kind=4), allocatable :: face_map(:)
  integer(kind=4), allocatable :: buf_face_map(:)
  integer                      :: cell_map_type
  integer                      :: bnd_cell_map_type
  integer                      :: face_map_type
  integer                      :: buf_face_map_type

  contains

  include 'Comm_Mod/Sequential/Allocate_Memory.f90'
  include 'Comm_Mod/Sequential/Close_File.f90'
  include 'Comm_Mod/Sequential/Create_New_Types.f90'
  include 'Comm_Mod/Sequential/End.f90'
  include 'Comm_Mod/Sequential/Exchange.f90'
  include 'Comm_Mod/Sequential/Global_Max_Real.f90'
  include 'Comm_Mod/Sequential/Global_Min_Real.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real.f90'
  include 'Comm_Mod/Sequential/Load_Buffers.f90'
  include 'Comm_Mod/Sequential/Load_Maps.f90'
  include 'Comm_Mod/Sequential/Open_File_Read.f90'
  include 'Comm_Mod/Sequential/Open_File_Write.f90'
  include 'Comm_Mod/Sequential/Read_Int.f90'
  include 'Comm_Mod/Sequential/Read_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Read_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Read_Face_Real.f90'
  include 'Comm_Mod/Sequential/Read_Text.f90'
  include 'Comm_Mod/Sequential/Start.f90'
  include 'Comm_Mod/Sequential/Wait.f90'
  include 'Comm_Mod/Sequential/Write_Int.f90'
  include 'Comm_Mod/Sequential/Write_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Write_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Write_Face_Real.f90'
  include 'Comm_Mod/Sequential/Write_Text.f90'

  end module
