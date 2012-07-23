function [ nodes, plainGrid, niceGrid ] = doPrecomputeGrid(x, psi, order)
  sdim = length(x);
  terms = length(psi);

  %
  % (*) The quadrature rule that is used here (Gauss-Hermite) assumes
  % the physicists' weight function. We need the probabilists' one,
  % therefore, a proper conversion is to be performed.
  %

  %
  % `level' means that for polynomial with of the order (2 * `level' - 1)
  % the integration will be exact. We want to have the exactness for
  % the order (2 * `order'), therefore, `level' = `order' + 1.
  %
  level = order + 1;

  [ nodes, weights, pointsSG ] = GaussianQuadrature.Physicists.constructSparseGrid(sdim, level);

  pointsTP = level^sdim;

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Type: Physicists' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', order }, ...
        { '  Accuracy level: %d', level }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  if pointsTP <= pointsSG
    [ nodes, weights ] = GaussianQuadrature.Physicists.constructTensorProduct(sdim, level);
  end

  points = min(pointsTP, pointsSG);

  %
  % See (*) above.
  %
  nodes = nodes * sqrt(2);

  plainGrid = zeros(terms, points);
  niceGrid = zeros(terms, points);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');

    plainGrid(i, :) = f(nodes);

    %
    % See (*) above.
    %
    niceGrid(i, :) = plainGrid(i, :) .* weights ./ pi^(sdim / 2);
  end
end
