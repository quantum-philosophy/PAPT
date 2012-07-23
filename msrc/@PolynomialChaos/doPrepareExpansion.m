function [ x, psi, index ] = doPrepareExpansion(sdim, order)
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
  %   * psi   - a vector of the Hermite polynomials,
  %   * index - the multi-index of the expansion.
  %

  terms = PolynomialChaos.countTerms(sdim, order);

  debug({ 'Construction of a new PC expansion.' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Number of terms: %d', terms });

  for i = 1:sdim
    x(i) = ipoly([ 'x', num2str(i) ]);
  end

  [ psi, index ] = PolynomialChaos.constructMD(x, order);

  assert(terms == length(psi), 'The number of terms is invalid.');
end
