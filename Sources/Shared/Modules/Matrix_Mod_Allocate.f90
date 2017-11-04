!==============================================================================!
  subroutine Matrix_Mod_Allocate(matrix, n_bnd_cells, n_cells, n_faces)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Matrix_Type) :: matrix
  integer           :: n_bnd_cells, n_cells, n_faces
!==============================================================================!

  ! Allocate memory for matrix
  allocate (matrix % row(n_cells+1));             matrix % row=0
  allocate (matrix % dia(n_cells));               matrix % dia=0
  allocate (matrix % sav(-n_bnd_cells:n_cells));  matrix % sav=0
  allocate (matrix % bou(-n_bnd_cells:-1));       matrix % bou=0
  allocate (matrix % pos(2,n_faces));             matrix % pos=0

  end subroutine
