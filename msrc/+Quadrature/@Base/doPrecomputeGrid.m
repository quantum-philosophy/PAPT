function [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(qd, x, psi, polynomialOrder, index)
  sdim = length(x);
  terms = length(psi);

  quadratureOrder = qd.polynomialOrderToQuadratureOrder(polynomialOrder);

  [ nodes, weights, pointsSG ] = qd.constructSparseGrid(sdim, quadratureOrder);

  pointsTP = qd.countTensorProductPoints(sdim, quadratureOrder);

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Type: Probabilists' }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', polynomialOrder }, ...
        { '  Quadrature order: %d', quadratureOrder }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  if pointsTP <= pointsSG
    [ nodes, weights ] = qd.constructTensorProduct(sdim, quadratureOrder);
  end

  points = min(pointsTP, pointsSG);

  plainGrid = zeros(terms, points);
  niceGrid = zeros(terms, points);

  norm = zeros(1, terms);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');
    plainGrid(i, :) = f(nodes);
    if nargin < 5
      norm(i) = sum(plainGrid(i, :).^2 .* weights);
    else
      norm(i) = prod(factorial(index(i, :) - 1));
    end
    niceGrid(i, :) = plainGrid(i, :) .* weights / norm(i);
  end

  assert(all(norm >= 0), 'Normalization constants cannot be negative.');
end
