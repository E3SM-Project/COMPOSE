! This file is modified from diag.F90 in
!     www.geosci-model-dev.net/5/887/2012/gmd-5-887-2012-supplment.zip
! supplementary material for the paper [LSPT2012]
!     Lauritzen, Skamarock, Prather, Taylor, A standard test case suite for
!     two-dimensional linear transport on the sphere, GMD 2012.

MODULE diag
  use iso_c_binding

  IMPLICIT NONE
  !
  ! Numerical mixing diagnostics are computed in Subroutine correlation_diag below.
  !
  ! Filament preservation diagnostic is computed in Subroutine filament_diag below.
  !
  ! It is assumed that the data is in one-dimensional arrays of length K, however,
  ! the code can easily be changed to accommodate other data structures.
  !
  ! dA(j) is the spherical area of grid cell j
  ! f1 is the mixing ratio of tracer 1
  ! f2 is the mixing ratio of tracer 2
  !

  integer, parameter :: dbl = c_double

CONTAINS
  !
  ! compute correlation diagnostics
  !

  SUBROUTINE correlation_diag(f1,f2,K,dA,real_mixing,overshooting,range_pres_unmixing) bind(c)
    IMPLICIT NONE
    INTEGER(c_int), value, INTENT(IN) :: K
    REAL(c_double), DIMENSION(K)      , INTENT(IN) :: f1,f2,dA
    real(c_double), intent(out) :: real_mixing,overshooting,range_pres_unmixing
    !
    ! local workspace
    !
    REAL(dbl)     :: root, tol,q1,q2,c
    INTEGER  :: j
    REAL(dbl)     :: q1_min,q1_max,q2_min,q2_max
    REAL(dbl)     :: total_area, sqrt_arg
    
    REAL(dbl), parameter :: eps = 1.0d-7

    q1_min = 0.1d0
    q1_max = 1.0d0

!! Changed so that end points of line segment are (q1_min,q2_min) and (q1_max,q2_max)
!! Note that "q2_min > q2_max".
    q2_min = corr_fct(q1_min)
    q2_max = corr_fct(q1_max)

    real_mixing          = 0.0d0
    overshooting         = 0.0d0
    range_pres_unmixing  = 0.0d0

    total_area = 0.0d0

    DO j = 1,K
       total_area = total_area+dA(j)

       q1 = f1(j)
       q2 = f2(j)

       !! Check to make sure we are not in the (very unlikely) situation where the argument
       !! to the sqrt (below) is negative.  This will happen if we are in a region which
       !! the cubic has three real roots, and so we'd have to be careful about which one we pick.
       sqrt_arg = (-DBLE(1687296) + DBLE(12168000)*q2 &
            -DBLE(29250000)*q2**2+DBLE(23437500)*q2**3+DBLE(29648025)*q1**2)
       IF (sqrt_arg < 0) THEN
          WRITE(6,*) 'Warning : (xk,yk) data is in region where there are three possible ', &
               ' real roots to the closest point problem'
          STOP
       ENDIF

       c = (DBLE(65340)*q1+12.d0*SQRT(-DBLE(1687296)+DBLE(12168000)*q2 &
            -DBLE(29250000)*q2**2+DBLE(23437500)*q2**3+DBLE(29648025)*q1**2))**(1.d0/3.d0)
       c=c/(DBLE(60))

       root = c-(-(DBLE(13)/DBLE(75))+(DBLE(5)/DBLE(12))*q2)/c
       root = MAX(0.1d0,root)
       root = MIN(1.0d0,root)

       !! Also fixed bug here : call to line_fct passed in q2 instead of q1 as argument.
       IF (q2 < corr_fct(q1)+eps.AND.q2 > line_fct(q1,q1_min,q1_max,q2_min,q2_max)-eps) THEN
          !
          ! `real' mixing
          !
          real_mixing = real_mixing + dist_fct(root,q1,q2)*dA(j)
       ELSE IF (q1 < q1_max+eps.AND.q1 > q1_min-eps.AND.q2 < q2_min+eps.AND.q2 > q2_max-eps) THEN
          !! Note that "q2_min > q2_max", so this 'if' branch had to be modifed

          !
          ! range-preserving unmixing
          !
          range_pres_unmixing = range_pres_unmixing+dist_fct(root,q1,q2)*dA(j)
       ELSE
          !
          ! overshooting
          !
          overshooting = overshooting + dist_fct(root,q1,q2)*dA(j)
       END IF
    END DO
    real_mixing = real_mixing/total_area
    range_pres_unmixing = range_pres_unmixing/total_area
    overshooting = overshooting/total_area
    return
    WRITE(*,*) "========================================================================"
    WRITE(*,*) " "
    WRITE(*,*) "mixing diagnostics"
    WRITE(*,*) " "
    WRITE(*,*) "------------------------------------------------------------------------"
    WRITE(*,*) " "
    WRITE(*,*) "real_mixing ",real_mixing
    WRITE(*,*) "range_pres_unmixing ",range_pres_unmixing
    WRITE(*,*) "overshooting     ",overshooting
    WRITE(*,*) " "
    WRITE(*,*) "========================================================================"
  END SUBROUTINE correlation_diag

  !
  ! correlation function
  !
  REAL(dbl) FUNCTION corr_fct(x)
    IMPLICIT NONE
    REAL(dbl) , INTENT(IN)  :: x
    corr_fct = -0.8d0*x**2+0.9d0
  END FUNCTION corr_fct
  !
  ! Eucledian distance function
  !
  REAL(dbl) FUNCTION dist_fct(x,x0,y0)
    IMPLICIT NONE
    REAL(dbl) , INTENT(IN)  :: x,x0,y0
    dist_fct = SQRT((x-x0)*(x-x0)/(0.9d0**2)+(corr_fct(x)-y0)*(corr_fct(x)-y0)/(0.792d0**2))
  END FUNCTION dist_fct
  !
  ! straight line line function
  !
  REAL(dbl) FUNCTION line_fct(x,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    REAL(dbl) , INTENT(IN)  :: x,xmin,xmax,ymin,ymax
    REAL(dbl)  :: a,b
    !
    ! line: y=a*x+b
    !
    a = (ymax-ymin)/(xmax-xmin)
    b = ymin-xmin*a
    line_fct = a*x+b
  END FUNCTION line_fct

  !
  ! linit = .TRUE. if t=0 else .FALSE.
  !
  ! Note that this subroutine must be called with the initial
  ! condition data to get "fila_t0".
  !
  SUBROUTINE filament_diag(K,f1,dA,fila_t0,linit,thresholds,fila_tf) bind(c)
    IMPLICIT NONE
    INTEGER(c_int), value, INTENT(IN)                    :: K
    REAL(c_double)   , DIMENSION(K)  , INTENT(IN)    :: f1, dA
    REAL(c_double)   , DIMENSION(100), INTENT(INOUT) :: fila_t0
    REAL(c_double)   , DIMENSION(100), INTENT(OUT)   :: thresholds, fila_tf
    INTEGER(c_int), value, INTENT(IN)    :: linit
    !
    ! local workspace
    !
    REAL(dbl)    :: threshold,tiny
    REAL(dbl)    :: out
    INTEGER :: j,jk,jlevels

    tiny=1.0d-12
    OPEN (unit = 31, file='filament.dat',status='replace')
    jlevels = 18
    IF (linit == 0) thresholds(jlevels+2) = -1
    DO jk=0,jlevels
       threshold = 0.1d0+(DBLE(jk)/DBLE(jlevels))*0.9d0
       IF (linit == 0) thresholds(jk+1) = threshold
       out       = 0.0d0
       DO j=1,K
          IF (f1(j).GE.threshold-tiny) THEN
             out = out+dA(j)
          END IF
       END DO
       IF (linit /= 0) THEN
          fila_t0(jk+1)=out
       ELSE
          IF (fila_t0(jk+1)<tiny) THEN
             fila_tf(jk+1) = 0.0d0
             !WRITE(31,*) threshold,0.0d0
          ELSE
             fila_tf(jk+1) = 100.0d0*out/fila_t0(jk+1)
             !WRITE(31,*) threshold,100.0d0*out/fila_t0(jk+1)
          END IF
       END IF
    END DO
    CLOSE(31)
  END SUBROUTINE filament_diag
END MODULE diag
