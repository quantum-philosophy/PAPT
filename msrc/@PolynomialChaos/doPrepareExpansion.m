function [ x, psi ] = doPrepareExpansion(sdim, order)
  %
  % Description:
  %
  %   Calculate multivariate (Hermite) polynomials of `sdim' r.v.'s
  %   of the maximal order `order'.
  %
  % Input:
  %
  %   * sdim   - the stochastic dimension of the expansion,
  %   * order  - the maximal order of the polynomials.
  %
  % Output:
  %
  %   * x     - a vector of `sdim' symbolic variables,
  %   * psi   - a vector of the Hermite polynomials.
  %

  terms = PolynomialChaos.countTerms(sdim, order);

  for i = 1:sdim
    x(i) = ipoly([ 'x', num2str(i) ]);
  end

  psi = PolynomialChaos.constructXD(x, terms);

  warn(sdim, order, terms);
end

function warn(sdim, order, terms)
  debug('------------------------------\n');
  debug('A new PC expansion is constructed:\n');
  debug('  Stochastic dimension: %d\n', sdim);
  debug('  Polynomial order: %d\n', order);
  debug('  Number of terms: %d\n', terms);
  debug('------------------------------\n');
end
