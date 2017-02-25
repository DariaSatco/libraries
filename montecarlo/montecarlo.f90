subroutine monte_carlo_nd ( func, dim_num, a, b, eval_num, result )

!*****************************************************************************80
!
!! MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
!
!  Discussion:
!
!    Unlike the other routines, this routine requires the user to specify
!    the number of function evaluations as an INPUT quantity.
!
!    No attempt at error estimation is made.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, a routine which evaluates
!    the function to be integrated, of the form:
!      function func ( x )
!      integer ( kind = 4 ) dim_num
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(dim_num)
!      func = ...
!      return
!      end
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the integration limits.
!
!    Input, integer ( kind = 4 ) EVAL_NUM, the number of function evaluations.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) eval_num
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(dim_num)

  result = 0.0D+00

  do i = 1, eval_num

    call random_vec( dim_num, a, b, x )

    result = result + func ( x )

  end do

  volume = product ( b(1:dim_num) - a(1:dim_num) )

  result = result * volume / real ( eval_num, kind = 8 )

  return
end

subroutine random_vec ( n, a, b, r )

!*****************************************************************************80
!
!! random_vec returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n) , a(n), b(n)

  do i = 1, n
    call random_number(r(i))
    r(i)=a(i)+(b(i)-a(i))*r(i)
  end do

  return
end
