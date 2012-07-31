function [ nodes, norm, grid, rvPower, rvProd, coeffMap ] = doPrepareExpansion(sdim, ddim, order, method)
  %
  % Input:
  %
  %   * sdim   - the stochastic dimension of the expansion,
  %   * sdim   - the deterministic dimension of the expansion,
  %   * order  - the maximal order of the polynomials.
  %

  terms = PolynomialChaos.countTerms(sdim, order);

  debug({ 'Construction of a new PC expansion.' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Number of terms: %d', terms });

  for i = 1:sdim
    x(i) = sympoly([ 'x', num2str(i) ]);
  end

  [ psi, index ] = PolynomialChaos.constructMD(x, order);

  assert(terms == length(psi), 'The number of terms is invalid.');

  %
  % Now, the Gaussian quadrature rule.
  %

  switch method
  case 'Physicists'
    gq = GaussianQuadrature.Physicists(x, psi, order);
  case 'Probabilists'
    gq = GaussianQuadrature.Probabilists(x, psi, order);
  case 'MultiProbabilists'
    gq = GaussianQuadrature.MultiProbabilists(x, psi, index);
  otherwise
    error('The method is unknown.');
  end

  points = gq.points;
  nodes = gq.nodes;
  norm = gq.norm;

  %
  % Precompute the grid for the given number of deterministic dimensions.
  %
  grid = zeros(ddim, points, terms);
  for i = 1:terms
    grid(:, :, i) = irep(gq.niceGrid(i, :), ddim, 1);
  end

  %
  % Construct the final polynomial with abstract coefficients
  % and produce the corresponding representation of it in terms
  % of monomials.
  %
  for i = 1:terms
    a(i) = sympoly([ 'a', num2str(i) ]);
  end

  [ rvPower, coeffMap ] = Utils.toMatrix(sum(a .* psi));

  mterms = size(rvPower, 2);

  assert(size(coeffMap, 1) == terms, 'The size of coeffMap is invalid.');
  assert(size(coeffMap, 2) == mterms, 'The size of coeffMap is invalid.');

  %
  % rvPower is a (sdim x mterms) matrix of the exponents of each of the monomials.
  %

  rvProd = zeros(mterms, points);

  for i = 1:mterms
    rvProd(i, :) = prod(realpow(nodes, irep(rvPower(:, i), 1, points)), 1);
  end
end
