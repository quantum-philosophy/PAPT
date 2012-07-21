function [ nodes, plainGrid, niceGrid ] = doPrecomputeGrid(x, psi, order)
  sdim = length(x);
  terms = length(psi);

  %
  % `level' means that for polynomial with of the order (2 * `level' - 1)
  % the integration will be exact. We want to have the exactness for
  % the order (2 * `order'), therefore, `level' = `order' + 1.
  %
  level = order + 1;

  [ nodes, weights, pointsSG ] = GaussianQuadrature.Peace.constructSparseGrid(sdim, level);

  pointsTP = level^sdim;

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Type: peace' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Accuracy level: %d', level }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  if pointsTP <= pointsSG
    [ nodes, weights ] = GaussianQuadrature.Peace.constructTensorProduct(sdim, level);
  end

  points = min(pointsTP, pointsSG);

  plainGrid = zeros(terms, points);
  niceGrid = zeros(terms, points);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');
    plainGrid(i, :) = f(nodes);
    niceGrid(i, :) = plainGrid(i, :) .* weights;
  end
end
