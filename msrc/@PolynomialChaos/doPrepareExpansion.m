function [ x, psi, norm, cq ] = doPrepareExpansion(sdim, order, level)
  %
  % Description:
  %
  %   Calculate multivariate (Hermite) polynomials of `sdim' r.v.'s
  %   of the maximal order `order'. In addition, compute the variances
  %   of each of the polynomials.
  %
  % Input:
  %
  %   * sdim   - the stochastic dimension of the expansion,
  %   * order  - the maximal order of the polynomials,
  %   * level  - the level of the Gauss quadrature to employ.
  %
  % Output:
  %
  %   * x     - a vector of `sdim' symbolic variables,
  %   * psi   - a vector of the Hermite polynomials,
  %   * norm  - a vector of the normalization coefficients of the expansion,
  %   * cq    - an instance of ChaosQuadrature to integrate with.
  %

  terms = PolynomialChaos.countTerms(sdim, order);

  for i = 1:sdim
    x(i) = ipoly([ 'x', num2str(i) ]);
  end

  psi = PolynomialChaos.constructXD(x, terms);

  cq = ChaosQuadrature(x, psi, order, level);

  norm = zeros(1, terms);
  norm(1) = 1;
  for i = 2:terms
    norm(i) = cq.integrateChaosProduct(i, i);
  end
end
