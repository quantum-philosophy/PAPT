function [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(x, psi, index)
  sdim = length(x);
  terms = length(psi);
  order = max(max(index));

  %
  % `level' means that for polynomial with of the order (2 * `level' - 1)
  % the integration will be exact. We want to have the exactness for
  % the order (2 * `order'), therefore, `level' = `order' + 1.
  %
  level = order + 1;

  [ nodes, weights, pointsSG ] = ...
    GaussianQuadrature.MultiProbabilists.constructSparseGrid(sdim, level);

  pointsTP = level^sdim;

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Type: MultiProbabilists' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Accuracy level: %d', level }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  if pointsTP <= pointsSG
    [ nodes, weights ] = ...
      GaussianQuadrature.MultiProbabilists.constructTensorProduct(sdim, level);
  end

  points = min(pointsTP, pointsSG);

  plainGrid = zeros(terms, points);
  niceGrid = zeros(terms, points);

  norm = zeros(1, terms);

  for i = 1:terms
    norm(i) = prod(factorial(index(i, :) - 1));
    f = Utils.toFunction(psi(i), x, 'rows');
    plainGrid(i, :) = f(nodes);
    niceGrid(i, :) = plainGrid(i, :) .* weights / norm(i);
  end
end
