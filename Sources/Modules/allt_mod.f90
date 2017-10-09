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

!----- Unknown type
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

!----- Matrix type
  type Matrix
    integer              :: nonzeros               ! number of nonzero entries
    real,    allocatable :: val(:)                 ! value
    real,    allocatable :: sav(:)                 ! saved value
    real,    allocatable :: bou(:)                 ! boundary value
    integer, allocatable :: row(:), col(:), dia(:) ! positions in the matrix
    integer, allocatable :: pos(:,:)               ! position in the matrix
  end type Matrix

end module 
