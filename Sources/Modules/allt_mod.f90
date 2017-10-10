!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!      Parameter definitions      !   Delft University of Technology   !
!       for all data types        !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module allt_mod

  implicit none

  !------------------!
  !   Unknown type   !
  !------------------!
  type Unknown          
    real, allocatable :: n(:)                ! new value
    real, allocatable :: o(:), oo(:)         ! old and older then old
    real, allocatable :: C(:), Co(:), Coo(:) ! convective fluxes
    real, allocatable :: Do(:), Doo(:)       ! difussive fluxes
    real, allocatable :: X(:), Xo(:), Xoo(:) ! surfce sources  
    real, allocatable :: mean(:)             ! long time average
    real, allocatable :: filt(:)             ! long time average
    real, allocatable :: q(:)                ! flux of a variable
    real, allocatable :: fluc(:) 
    real              :: URF                 ! under relaxation factor
    real              :: Stol                ! solver tolerance
    real              :: bound(128)          ! boundary values
    real              :: init(128)           ! initial values
    real              :: pro(11024)          ! inlfow profile
    real              :: Sigma               ! sigma
  end type Unknown

  !--------------------------------------------------------------------!
  !   Matrix type                                                      !
  !                                                                    !
  !   Matrix is stored in compressed row format.                       !
  !   (See: http://netlib.org/linalg/html_templates/node91.html)       !
  !                                                                    !
  !   Example:                                                         !
  !                                                                    !
  !       c   c  .    c                                                !
  !       o   o  .    o                                                !
  !       l   l       l                                                !
  !                                                                    !
  !       1   2       n                                                !
  !                                                                    !
  !    [ 10   0   4   5 ]  --> row 1                                   !
  !    [  2  12  -1   0 ]  --> rows store discretized control volumes  !
  !    [  0   1  99   7 ]  ...                                         !
  !    [ -3  11   0  53 ]  --> row n                                   !
  !                                                                    !
  !   Compressed row storage of the above matrix reads:                !
  !                                                                    !
  !   A % val = [  10   4   5   2  12  -1   1  99   7  -3  11  53 ]    !
  !   A % col = [   1   3   4   1   2   3   2   3   4   1   2   4 ]    !
  !   A % row = [   1   4   7  10 ]                                    !
  !                                                                    !
  !   A % dia = [   1   5   9  12 ]                                    !
  !--------------------------------------------------------------------!
  type Matrix
    integer              :: nonzeros               ! number of nonzero entries
    real,    allocatable :: val(:)                 ! value
    real,    allocatable :: sav(:)                 ! saved value
    real,    allocatable :: bou(:)                 ! boundary value
    integer, allocatable :: row(:), col(:), dia(:) ! positions in the matrix
    integer, allocatable :: pos(:,:)               ! position in the matrix
  end type Matrix

end module 
