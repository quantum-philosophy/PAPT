function sandia_rules_test35 ( )

%*****************************************************************************80
%
%% SANDIA_RULES_TEST35 uses HERMITE_INTEPOLAN using Chebyshev spaced data.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 November 2011
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SANDIA_RULES_TEST35\n' );
  fprintf ( 1, '  HERMITE_INTERPOLANT computes the Hermite interpolant to data.\n' );
  fprintf ( 1, '  Here, f(x) is the Runge function\n' );
  fprintf ( 1, '  and the data is evaluated at Chebyshev spaced points.\n' );
  fprintf ( 1, '  As N increases, the maximum error decreases.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '     N     Max | F(X) - H(F(X)) |\n' );
  fprintf ( 1, '\n' );

  for n = 3 : 2 : 15

    nd = 2 * n;
    ns = 10 * ( n - 1 ) + 1;
    xlo = -5.0;
    xhi = +5.0;
    x = r8vec_chebyshev ( n, xlo, xhi );

    y(1:n) = 1.0 ./ ( 1.0 + x(1:n).^2 );
    yp(1:n) = - 2.0 * x(1:n) ./ ( 1.0 + x(1:n).^2 ).^2;

    [ xd, yd, xdp, ydp ] = hermite_interpolant ( n, x, y, yp );
%
%  Compare exact and interpolant at sample points.
%
    xs = linspace ( xlo, xhi, ns );

    [ ys, ysp ] = hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs );

    max_dif = 0.0;
    for i = 1 : ns
      xt = xs(i);
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = max ( max_dif, abs ( ys(i) - yt ) );
    end

    fprintf ( 1, '  %4d  %14g\n', n, max_dif );

  end

  return
end
