function [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(x, psi, index)
  sdim = length(x);
  terms = length(psi);
  order = PolynomialChaos.indexToOrder(index);

  %
  % `order' of a Gaussian quadrature rule means that for polynomials with
  % order (2 * `order' - 1) the integration will be exact. We want to have
  % the exactness for the order (2 * `order'), therefore, `order' := `order' + 1.
  %
  requiredOrder = order + 1;

  [ nodes, weights, pointsSG ] = ...
    GaussianQuadrature.MultiProbabilists.constructSparseGrid(sdim, requiredOrder);

  pointsTP = requiredOrder^sdim;

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Type: MultiProbabilists' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  if pointsTP <= pointsSG
    [ nodes, weights ] = ...
      GaussianQuadrature.MultiProbabilists.constructTensorProduct(sdim, requiredOrder);
  end

  points = min(pointsTP, pointsSG);

  plainGrid = zeros(terms, points);
  niceGrid = zeros(terms, points);

  norm = zeros(1, terms);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');
    plainGrid(i, :) = f(nodes);
    norm(i) = prod(factorial(index(i, :) - 1));
    niceGrid(i, :) = plainGrid(i, :) .* weights / norm(i);
  end
end
