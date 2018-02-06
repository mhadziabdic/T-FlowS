!==============================================================================!
      subroutine Sort_Int_Carry_Int(x, y, n, kflag)
!------------------------------------------------------------------------------!
!   Sorts int. array x and makes the same changes in int. array y.             !
!   It was downladed from nist and then slightly modified.                     !
!------------------------------------------------------------------------------!
      dimension il(21), iu(21)
      integer   x(n), y(n), t, tt, ty, tty
!==============================================================================!

!
! First executable statement 
!
      nn = n
      if (nn.ge.1) go to 10
      write(*,*) 'isort: the number of values to be sorted was not positive.'
      return
   10 kk = iabs(kflag)
      if ((kk.eq.1).or.(kk.eq.2)) go to 15
      write(*,*) 'isort: the sort control parameter, k, was not 2, 1, -1, or -2.'
      return
!
! Alter array x to get decreasing order if needed
!
   15 if (kflag.ge.1) go to 30
      do 20 i=1,nn
      x(i) = -x(i)
   20 continue
   30 if(kk==1) goto 100
      if(kk==2) goto 200 
!
! Sort x only
!
  100 continue
      m=1
      i=1
      j=nn
      r=.375
  110 if (i .eq. j) go to 155
      if (r .gt. .5898437) go to 120
      r=r+3.90625e-2
      go to 125
  120 r=r-.21875
  125 k=i
!                                  select a central element of the
!                                  array and save it in location t
      ij = i + ifix (float (j-i) * r)
      t=x(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) go to 130
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
  130 l=j
!                                  if last element of array is less than
!                                  t, interchange with t
      if (x(j) .ge. t) go to 140
      x(ij)=x(j)
      x(j)=t
      t=x(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) go to 140
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
      go to 140
  135 tt=x(l)
      x(l)=x(k)
      x(k)=tt
!                                  find an element in the second half of
!                                  the array which is smaller than t
  140 l=l-1
      if (x(l) .gt. t) go to 140
!                                  find an element in the first half of
!                                  the array which is greater than t
  145 k=k+1
      if (x(k) .lt. t) go to 145
!                                  interchange these elements
      if (k .le. l) go to 135
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      if (l-i .le. j-k) go to 150
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 160
  150 il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 160
!                                  begin again on another portion of
!                                  the unsorted array
  155 m=m-1
      if (m .eq. 0) go to 300
      i=il(m)
      j=iu(m)
  160 if (j-i .ge. 1) go to 125
      if (i .eq. 1) go to 110
      i=i-1
  165 i=i+1
      if (i .eq. j) go to 155
      t=x(i+1)
      if (x(i) .le. t) go to 165
      k=i
  170 x(k+1)=x(k)
      k=k-1
      if (t .lt. x(k)) go to 170
      x(k+1)=t
      go to 165
!
! Sort x and carry y along
!
  200 continue
      m=1
      i=1
      j=nn
      r=.375
  210 if (i .eq. j) go to 255
      if (r .gt. .5898437) go to 220
      r=r+3.90625e-2
      go to 225
  220 r=r-.21875
  225 k=i
!                                  select a central element of the
!                                  array and save it in location t
      ij = i + ifix (float (j-i) *r)
      t=x(ij)
      ty= y(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) go to 230
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
       y(ij)= y(i)
       y(i)=ty
      ty= y(ij)
  230 l=j
!                                  if last element of array is less than
!                                  t, interchange with t
      if (x(j) .ge. t) go to 240
      x(ij)=x(j)
      x(j)=t
      t=x(ij)
       y(ij)= y(j)
       y(j)=ty
      ty= y(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) go to 240
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
       y(ij)= y(i)
       y(i)=ty
      ty= y(ij)
      go to 240
  235 tt=x(l)
      x(l)=x(k)
      x(k)=tt
      tty= y(l)
       y(l)= y(k)
       y(k)=tty
!                                  find an element in the second half of
!                                  the array which is smaller than t
  240 l=l-1
      if (x(l) .gt. t) go to 240
!                                  find an element in the first half of
!                                  the array which is greater than t
  245 k=k+1
      if (x(k) .lt. t) go to 245
!                                  interchange these elements
      if (k .le. l) go to 235
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      if (l-i .le. j-k) go to 250
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 260
  250 il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 260
!                                  begin again on another portion of
!                                  the unsorted array
  255 m=m-1
      if (m .eq. 0) go to 300
      i=il(m)
      j=iu(m)
  260 if (j-i .ge. 1) go to 225
      if (i .eq. 1) go to 210
      i=i-1
  265 i=i+1
      if (i .eq. j) go to 255
      t=x(i+1)
      ty= y(i+1)
      if (x(i) .le. t) go to 265
      k=i
  270 x(k+1)=x(k)
       y(k+1)= y(k)
      k=k-1
      if (t .lt. x(k)) go to 270
      x(k+1)=t
       y(k+1)=ty
      go to 265
!
! Clean up
!
  300 if (kflag.ge.1) return
      do 310 i=1,nn
      x(i) = -x(i)
  310 continue
      return
      end subroutine
