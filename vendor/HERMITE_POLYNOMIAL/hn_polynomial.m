function p = hn_polynomial ( m, n, x )

%*****************************************************************************80
%
%% HN_POLYNOMIAL evaluates the normalized physicist's Hermite polynomials.
%
%  Discussion:
%
%    These polynomials satisfy the orthonormality condition:
%
%      Integral ( -oo < X < +oo ) exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    25 February 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz, Irene Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964,
%    ISBN: 0-486-61272-4,
%    LC: QA47.A34.
%
%    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
%    NIST Handbook of Mathematical Functions,
%    Cambridge University Press, 2010,
%    ISBN: 978-0521192255,
%    LC: QA331.N57.
%
%  Parameters:
%
%    Input, integer M, the number of evaluation points.
%
%    Input, integer N, the highest order polynomial to compute.
%    Note that polynomials 0 through N will be computed.
%
%    Input, real X(M,1), the evaluation points.
%
%    Output, real P(M,N+1), the values of the polynomials of index 0 through N.
%
  p = zeros ( m, n + 1 );

  p(1:m,1) = 1.0;

  if ( n == 0 )
    return
  end

  p(1:m,2) = 2.0 * x(1:m,1);
 
  for j = 2 : n
    p(1:m,j+1) = 2.0 * x(1:m,1) .* p(1:m,j) - 2.0 * ( j - 1 ) * p(1:m,j-1);
  end
%
%  Normalize.
%
  d = diag ( 1.0 ./ sqrt ( gamma ( 1:n+1) .* 2.^(0:n) * sqrt ( pi ) ) );

  p = p * d;

  return
end
