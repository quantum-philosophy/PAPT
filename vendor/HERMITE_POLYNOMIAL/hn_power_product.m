function table = hn_power_product ( p, e )

%*****************************************************************************80
%
%% HN_POWER_PRODUCT: normalized physicist's Hermite power products.
%
%  Discussion:
%
%    Let Hn(I,X) represent the normalized physicist's Hermite  polynomial 
%    of degree I.  
%
%    For polynomial chaos applications, it is of interest to know the
%    value of the integrals of products of X with every possible pair
%    of basis functions.  That is, we'd like to form
%
%      Tij = Integral ( -oo < X < +oo ) X^E * Hn(I,X) * Hn(J,X) exp(-X*X) dx
%
%    We will estimate these integrals using Gauss-Hermite quadrature.
%
%    When E is 0, the computed table should be the identity matrix.
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
%    Input, integer P, the maximum degree of the polyonomial factors.
%    0 <= P.
%
%    Input, integer E, the exponent of X in the integrand.
%    0 <= E.
%
%    Output, real TABLE(P+1,P+1), the table of integrals.  TABLE(I,J)
%    represents the weighted integral of X^E * Hn(I+1,X) * Hn(J+1,X).
%
  table(1:p+1,1:p+1) = 0.0;

  order = p + 1 + floor ( ( e + 1 ) / 2 );

  [ x_table, w_table ] = h_quadrature_rule ( order );

  for k = 1 : order

    x = x_table(k);
    h_table = hn_polynomial ( 1, p, x );
%
%  The following formula is an outer product in H_TABLE.
%
    if ( e == 0 )
      table(1:p+1,1:p+1) = table(1:p+1,1:p+1) ...
        + w_table(k) * h_table(1:p+1)' * h_table(1:p+1);
    else
      table(1:p+1,1:p+1) = table(1:p+1,1:p+1) ...
        + w_table(k) * x^e * h_table(1:p+1)' * h_table(1:p+1);
    end

  end

  return
end

