function hermite_polynomial_test10 ( p, b )

%*****************************************************************************80
%
%% HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.
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
%  Parameters:
%
%    Input, integer P, the maximum degree of the polynomial factors.
%
%    Input, real B, the coefficient of X in the exponential factor.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'HERMITE_POLYNOMIAL_TEST10\n' );
  fprintf ( 1, '  Compute a normalized probabilist''s Hermite exponential product table.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-0.5*X*X) dx\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  where Hen(I,X) = normalized physicist''s Hermite polynomial of degree I.\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Maximum degree P = %d\n', p );
  fprintf ( 1, '  Exponential argument coefficient B = %f\n', b );

  table = hen_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' );

  return
end
