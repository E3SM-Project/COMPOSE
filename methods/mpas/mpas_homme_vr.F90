! This file is derived from
!     https://github.com/E3SM-Project/E3SM/blob/master/components/homme/src/share/vertremap_base.F90
! which is licensed according to
!     https://github.com/E3SM-Project/E3SM/blob/master/LICENSE

module vertremap_base

  implicit none

  integer, parameter :: real_kind = 8, r8 = real_kind, int_kind = 4

contains

  subroutine remap_Q_ppm(Qdp,nx,qsize,nlev,vert_remap_q_alg,dp1,dp2)
    ! remap 1 field
    ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
    !         dp1   layer thickness (source)
    !         dp2   layer thickness (target)
    ! vert_remap_q_alg
    !   1  mirrored values in ghost cells (i.e., no flux)
    !   2  piecewise constant in 2 cells near boundaries (don't use ghost cells)
    !   3  UKMO splines ghost cells (copy last value into all ghost cells)
    !  10  linear extrap with limiter
    !  11  linear extrap
    !
    ! output: remaped Qdp, conserving mass
    !
    implicit none
    integer,intent(in) :: nx,qsize,nlev,vert_remap_q_alg
    real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
    real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
    ! Local Variables
    integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
    real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
    real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
    real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
    real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
    real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
    real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
    real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
    real(kind=real_kind), dimension(       nlev   ) :: z1, z2
    real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
    real(kind=real_kind) :: mymass, massn1, massn2, ext(2)
    integer :: i, j, k, q, kk, kid(nlev)

    if (vert_remap_q_alg < 1 .or. vert_remap_q_alg > 11 .or. &
         (vert_remap_q_alg > 3 .and. vert_remap_q_alg < 9)) then
       print *, 'remap_Q_ppm INPUT ERROR: alg', vert_remap_q_alg, 'is invalid'
       return
    end if

    do j = 1 , nx
       do i = 1 , nx

          pin(1)=0
          pio(1)=0
          do k=1,nlev
             dpn(k)=dp2(i,j,k)
             dpo(k)=dp1(i,j,k)
             pin(k+1)=pin(k)+dpn(k)
             pio(k+1)=pio(k)+dpo(k)
          enddo

          pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
          !It makes sure there's an old interface value below the domain that is larger.
          pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
          !Therefore, the pressure of that mass cannot either.
          !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
          do k = 1 , gs
             dpo(1   -k) = dpo(       k)
             dpo(nlev+k) = dpo(nlev+1-k)
          enddo

          !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
          !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
          !interface location. In practice, the grid never deforms past one cell, so the search can be
          !simplified by this. Also, the interval of integration is usually of magnitude close to zero
          !or close to dpo because of minimial deformation.
          !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
          !I set them equal to each other.
          do k = 1, nlev
             kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
             !Find the index of the old grid cell in which this new cell's bottom interface resides.
             if (pio(kk) <= pin(k+1)) then
                do while ( pio(kk) <= pin(k+1) )
                   kk = kk + 1
                enddo
                kk = kk - 1                   !kk is now the cell index we're integrating over.
             else
                call binary_search(pio, pin(k+1), kk, nlev)
             end if
             if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
             !Top bounds match anyway, so doesn't matter what coefficients are used
             kid(k) = kk                   !Save for reuse
             z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
             !In fact, we're usually integrating very little or almost all of the cell in question
             z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
             !coordinate domain [-0.5,0.5].
          enddo

          !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
          !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
          ppmdx(:,:) = compute_ppm_grids( dpo, nlev )

          !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
          !mass accumulation
          do q = 1 , qsize
             !Accumulate the old mass up to old grid cell interface locations to simplify integration
             !during remapping. Also, divide out the grid spacing so we're working with actual tracer
             !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
             !tracer consistency for an initially uniform field. I copied it from the old remap routine.
             masso(1) = 0.
             do k = 1 , nlev
                ao(k) = Qdp(i,j,k,q)
                masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
                ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
             enddo
             if (vert_remap_q_alg >= 10) then
                ext(1) = minval(ao(1:nlev))
                ext(2) = maxval(ao(1:nlev))
                call linextrap(dpo(2), dpo(1), dpo(0), dpo(-1), &
                     ao(2), ao(1), ao(0), ao(-1), &
                     vert_remap_q_alg, ext(1), ext(2))
                call linextrap(dpo(nlev-1), dpo(nlev), dpo(nlev+1), dpo(nlev+2), &
                     ao(nlev-1), ao(nlev), ao(nlev+1), ao(nlev+2), &
                     vert_remap_q_alg, ext(1), ext(2))
             else
                !Fill in ghost values. Ignored if vert_remap_q_alg == 2
                do k = 1 , gs
                   if (vert_remap_q_alg == 3) then
                      ao(1   -k) = ao(1)
                      ao(nlev+k) = ao(nlev)
                   elseif (vert_remap_q_alg == 2 .or. vert_remap_q_alg == 1) then   !Ignored if vert_remap_q_alg == 2
                      ao(1   -k) = ao(       k)
                      ao(nlev+k) = ao(nlev+1-k)
                   endif
                enddo
             end if
             !Compute monotonic and conservative PPM reconstruction over every cell
             coefs(:,:) = compute_ppm( ao , ppmdx, nlev, vert_remap_q_alg )
             !Compute tracer values on the new grid by integrating from the old cell bottom to the new
             !cell interface to form a new grid mass accumulation. Taking the difference between
             !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
             !supposed to hold the full mass this needs no normalization.
             massn1 = 0.
             do k = 1 , nlev
                kk = kid(k)
                massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
                Qdp(i,j,k,q) = massn2 - massn1
                massn1 = massn2
             enddo
          enddo
       enddo
    enddo
  end subroutine remap_Q_ppm

  !THis compute grid-based coefficients from Collela & Woodward 1984.
  function compute_ppm_grids( dx, nlev )   result(rslt)
    implicit none
    real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
    integer, intent(in) :: nlev
    real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
    integer :: j
    integer :: indB, indE

    !Calculate grid-based coefficients for stage 1 of compute_ppm
    do j = 0 , nlev+1
       rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
       rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
       rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
    enddo

    !Caculate grid-based coefficients for stage 2 of compute_ppm
    do j = 0 , nlev
       rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
       rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
       rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
       rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
       rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
       rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
       rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
    enddo
  end function compute_ppm_grids

  !This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
  function compute_ppm( a , dx, nlev, vert_remap_q_alg )    result(coefs)
    implicit none
    real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
    real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
    integer, intent(in) :: nlev, vert_remap_q_alg
    real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
    real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
    real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
    real(kind=real_kind) :: da                                !Ditto
    ! Hold expressions based on the grid (which are cumbersome).
    real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
    real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
    integer :: j
    integer :: indB, indE

    ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
    do j = 0 , nlev+1
       da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
       dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1._real_kind,da)
       if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
    enddo

    ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
    do j = 0 , nlev
       ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
            * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
    enddo

    ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
    ! (dimension nlev) using ai on either side and ao within the cell.
    do j = 1 , nlev
       al = ai(j-1)
       ar = ai(j  )
       if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
          al = a(j)
          ar = a(j)
       endif
       if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
       if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
       !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
       coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
       coefs(1,j) = ar - al
       ! coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
       coefs(2,j) = 3. * (-2. * a(j) + ( al + ar ))
    enddo

    !If vert_remap_q_alg == 2, use piecewise constant in the boundaries, and don't use ghost cells.
    if (vert_remap_q_alg == 2) then
       coefs(0,1:2) = a(1:2)
       coefs(1:2,1:2) = 0.
       coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
       coefs(1:2,nlev-1:nlev) = 0.D0
    endif
  end function compute_ppm

  !Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
  !given two bounds. Make sure this gets inlined during compilation.
  function integrate_parabola( a , x1 , x2 )    result(mass)
    implicit none
    real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
    real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
    real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
    real(kind=real_kind)             :: mass
    mass = a(0) * (x2 - x1) + a(1) * (x2 * x2 - x1 * x1) / 0.2D1 + a(2) * (x2 * x2 * x2 - x1 * x1 * x1) / 0.3D1
  end function integrate_parabola

  ! Find k such that pio(k) <= pivot < pio(k+1). Provide a reasonable input
  ! value for k.
  subroutine binary_search(pio, pivot, k, nlev)
    real(kind=real_kind), intent(in) :: pio(nlev+2), pivot    
    integer, intent(inout) :: k
    integer, intent(in) :: nlev
    integer :: lo, hi, mid

    if (pio(k) > pivot) then
       lo = 1
       hi = k
    else
       lo = k
       hi = nlev+2
    end if
    do while (hi > lo + 1)
       k = (lo + hi)/2
       if (pio(k) > pivot) then
          hi = k
       else
          lo = k
       end if
    end do
    k = lo
  end subroutine binary_search

  subroutine linextrap(dx1,dx2,dx3,dx4,y1,y2,y3,y4,q_alg,lo,hi)
    integer, intent(in) :: q_alg
    real(kind=real_kind), intent(in) :: dx1,dx2,dx3,dx4,y1,y2,lo,hi
    real(kind=real_kind), intent(out) :: y3,y4

    real(kind=real_kind) :: den,num,a
    real(kind=real_kind) :: z3,z4

    ! In exact arithmetic, the following is equivalent to
    !   x1 = half*dx1
    !   x2 = x1 + half*(dx1 + dx2)
    !   x3 = x2 + half*(dx2 + dx3)
    !   x4 = x3 + half*(dx3 + dx4)
    !   a  = (x3-x1)/(x2-x1)
    !   y3 = (1-a)*y1 + a*y2
    !   a  = (x4-x1)/(x2-x1)
    !   y4 = (1-a)*y1 + a*y2
    ! In FP, -fp-model-fast -O3 produces different results depending
    ! on whether -qopenmp is present.

    den = (dx1 + dx2)/2
    num = den + (dx2 + dx3)/2
    a = num/den
    y3 = (1-a)*y1 + a*y2

    num = num + (dx3 + dx4)/2
    a  = num/den
    y4 = (1-a)*y1 + a*y2

    if (q_alg /= 10) return
    y3 = max(lo, min(hi, y3))
    y4 = max(lo, min(hi, y4))
  end subroutine linextrap

  subroutine eval_lagrange_poly(n, xs, ys, xi, y)
    integer, intent(in) :: n
    real(kind=real_kind), intent(in) :: xs(:), ys(:), xi
    real(kind=real_kind), intent(out) :: y

    integer :: i, j
    real(kind=real_kind) :: f

    y = 0
    do i = 1,n
       f = 1
       do j = 1,n
          if (i == j) cycle
          f = f*((xi - xs(j))/(xs(i) - xs(j)))
       end do
       y = y + ys(i)*f
    end do
  end subroutine eval_lagrange_poly

  ! ----------------------------------------------------------------------------
  ! Testing

  function eval_f(x) result(f)
    real(real_kind), intent(in) :: x
    real(real_kind) :: f

    f = cos(7.d0*x**2 + 0.33) + x**3 - 0.4*x**2 + x
    !                                           ^ this term makes the OOA 1.5, except for linextrap
  end function eval_f

  subroutine test()
    real(real_kind), parameter :: pi = 3.141592653589793d0
    integer, parameter :: nme = 12, nlev_max = 2**nme
    integer, parameter :: algs(6) = (/1,2,3,10,11,9/)

    real(real_kind) :: p1(nlev_max+1), p2(nlev_max+1), dp1(nlev_max), dp2(nlev_max), &
         Qdp(1,1,nlev_max,1), c1, re, re_prev, ooa, f, ft, num, den
    integer :: alg, algidx, i, k, nlev
    logical :: first

    c1 = 0.5d0

    do algidx = 1, size(algs)
       alg = algs(algidx)
       write(*,'(a,i3)') 'alg>',alg
       first = .true.
       do i = 3,nme
          nlev = 2**i

          p1(1) = 0
          do k = 2,nlev+1
             p1(k) = sin(c1*pi*real(k,real_kind)/(nlev+1))/sin(c1*pi)
             dp1(k-1) = p1(k) - p1(k-1)
          end do

          p2(1) = 0
          do k = 2,nlev+1
             dp2(k-1) = 1.0d0/nlev
             p2(k) = p2(k-1) + dp2(k-1)
          end do

          do k = 1,nlev
             Qdp(1,1,k,1) = dp1(k)*eval_f(p1(k) + 0.5_r8*dp1(k))
          end do

          call remap_Q_ppm(Qdp, 1, 1, nlev, alg, dp1, dp2)

          ! Point evaluation rather than integrating over the cell
          ! limits OOA to 2.
          num = 0; den = 0
          do k = 1,nlev
             ft = eval_f(p2(k) + 0.5_r8*dp2(k))
             f = Qdp(1,1,k,1)/dp2(k)
             num = num + abs(f - ft)**2
             den = den + ft**2
          end do
          re = sqrt(num/den)

          if (first) then
             ooa = -1
          else
             ooa = -(log(re) - log(re_prev))/log(2.0_r8)
          end if
          write(*,'(a,i2,es12.2,f6.2)') 're> ',i,re,ooa
          re_prev = re
          first = .false.
       end do
    end do
  end subroutine test

  subroutine remap1c(nlev, alg, dp1i, dp2i, Qdpi) bind(c)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: nlev, alg
    real(kind=c_double), intent(in) :: dp1i(nlev), dp2i(nlev)
    real(kind=c_double), intent(inout) :: Qdpi(nlev)

    real(kind=real_kind) :: dp1(1,1,nlev), dp2(1,1,nlev), Qdp(1,1,nlev,1)

    dp1(1,1,:nlev) = dp1i(:nlev)
    dp2(1,1,:nlev) = dp2i(:nlev)
    Qdp(1,1,:nlev,1) = Qdpi(:nlev)
    call remap_Q_ppm(Qdp, 1, 1, nlev, alg, dp1, dp2)
    Qdpi(:nlev) = Qdp(1,1,:nlev,1)
  end subroutine remap1c

end module vertremap_base

#ifdef VR_MAIN
program main
  ! standalone conv tester:
  !   gfortran -DVR_MAIN mpas_homme_vr.F90; if [ $? == 0 ]; then ./a.out; fi
  ! shared lib for use in python:
  !   gfortran -c -fPIC mpas_homme_vr.F90; gfortran -shared mpas_homme_vr.o -o libvr.so

  use vertremap_base, only: test

  call test()
end program main
#endif
