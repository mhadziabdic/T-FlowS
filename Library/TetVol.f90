!======================================================================!
  REAL FUNCTION TetVol(xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD)
!----------------------------------------------------------------------!
!   Returns the volume of tethraedra spanned with A, B, C and D.       !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  REAL :: xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD
!======================================================================!
!                                                                      !
!   The order of nodes (A,B,C and D) DOES matters.                     !
!                                                                      !
!                D-----C                                               !
!               / \  . |                                               !
!              /   \   |                                               !
!             /  .  \  |    I am not 100% sure that the figure is OK   !
!            / .     \ |                                               !
!           /.        \|                                               !
!          A-----------B                                               !
!                                                                      !
!----------------------------------------------------------------------!

  TetVol=( ((yB-yA)*(zC-zA)-(yC-yA)*(zB-zA))*(xD-xA) +              &
           ((xC-xA)*(zB-zA)-(xB-xA)*(zC-zA))*(yD-yA) +              &
           ((xB-xA)*(yC-yA)-(xC-xA)*(yB-yA))*(zD-zA) )/6.0

  END FUNCTION TetVol
