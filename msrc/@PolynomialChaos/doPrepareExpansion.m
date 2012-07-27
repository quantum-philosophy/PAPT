function [ nodes, norm, grid, coeffMap, rvProd ] = doPrepareExpansion(sdim, ddim, order)
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

  gq = GaussianQuadrature.MultiProbabilists(x, psi, index);

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

  [ P, coeffMap ] = Utils.toMatrix(sum(a .* psi));

  %
  % P is a (sdim x mterms) matrix of the exponents of each of the monomials.
  %

  mterms = size(P, 2);
  rvProd = zeros(mterms, points);

  for i = 1:mterms
    rvProd(i, :) = prod(realpow(nodes, irep(P(:, i), 1, points)), 1);
  end
end
