function [ nodes, norm, grid, rvPower, rvProd, coeffMap ] = ...
  doPrepareExpansion(pc, sdim, ddim, order, method)
  %
  % Input:
  %
  %   * sdim   - the stochastic dimension of the expansion,
  %   * sdim   - the deterministic dimension of the expansion,
  %   * order  - the maximal order of the polynomials.
  %

  terms = pc.countTerms(sdim, order);

  debug({ 'Construction of a new PC expansion.' }, ...
        { '  Method: %s', PolynomialChaos.methodName(method) }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Number of terms: %d', terms });

  for i = 1:sdim
    x(i) = sympoly([ 'x', num2str(i) ]);
  end

  [ psi, index ] = pc.constructMD(x, order);

  assert(terms == length(psi), 'The number of terms is invalid.');

  %
  % Now, the Gaussian quadrature rule.
  %

  qd = Quadrature.(method.quadratureName)(x, psi, order, index, method);

  points = qd.points;
  nodes = qd.nodes;
  norm = qd.norm;

  %
  % Precompute the grid for the given number of deterministic dimensions.
  %
  grid = cell(terms);
  for i = 1:terms
    grid{i} = irep(qd.niceGrid(i, :), ddim, 1);
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
