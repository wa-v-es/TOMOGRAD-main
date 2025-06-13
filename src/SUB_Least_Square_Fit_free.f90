SUBROUTINE LINEARFIT(x,y,npts,slope_best,ypoint,std)
! Get the linear fit for the given x y series and spit out slope  and intersection
! Assuming the first point is always at x=0 and the line must go through it
!======================================================================================
IMPLICIT NONE
INTEGER,PARAMETER :: maxp = 1000000
REAL, DIMENSION(maxp) :: x,y,slope_initial,ypoint_initial,rms,slope
INTEGER :: npts,num_step,i,j,k,l,m,n
REAL :: ypoint,slope_step,slope_best,value,max_slope,min_slope,ave_slope,std
REAL :: a1,a2,b0,b1,d,a,b,d1
!======================================================================================

  a1 = 0
  a2 = 0
  b0 = 0
  b1 = 0
  DO i = 1, npts
    a1 = a1 + x(i)
    a2 = a2 + x(i) * x(i)
    b0 = b0 + y(i)
    b1 = b1 + y(i) * x(i)
  END DO
  a1 = a1 / npts
  a2 = a2 / npts
  b0 = b0 / npts
  b1 = b1 / npts
  d = a1 * a1 - a2
  a = a1 * b1 - a2 * b0
  a = a / d
  b = a1 * b0 - b1
  b = b / d
!  Evaluation of standard deviation d (unbiased estimate) 
  d = 0
  DO i = 1, npts
    d1 = y(i) - a - b * x(i)
    d  = d + d1 * d1
  END DO
  d = SQRT(d / npts)
slope_best=b
ypoint=a
std=d
END SUBROUTINE LINEARFIT
