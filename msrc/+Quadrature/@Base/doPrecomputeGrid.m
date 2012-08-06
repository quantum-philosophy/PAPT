function [ nodes, grid, norm ] = doPrecomputeGrid(qd, x, psi, polynomialOrder, index, method)

  sdim = length(x);
  terms = length(psi);

  quadratureOrder = qd.polynomialOrderToQuadratureOrder(polynomialOrder);

  pointsSG = qd.countSparseGridPoints(sdim, quadratureOrder);
  pointsTP = qd.countTensorProductPoints(sdim, quadratureOrder);

  debug({ 'Precomputation of a new grid.' }, ...
        { '  Method: %s', Quadrature.methodName(method) }, ...
        { '  Stochastic dimensions: %d', sdim }, ...
        { '  Polynomial order: %d', polynomialOrder }, ...
        { '  Quadrature order: %d', quadratureOrder }, ...
        { '  Number of terms: %d', terms }, ...
        { '  Sparse grid points: %d', pointsSG }, ...
        { '  Tensor product points: %d', pointsTP });

  type = lower(method.quadratureType);

  if strcmp(type, 'adaptive')
    if pointsTP <= pointsSG
      type = 'tensor';
    else
      type = 'sparse';
    end
  end

  switch type
  case 'sparse'
    [ nodes, weights, points ] = qd.constructSparseGrid(sdim, quadratureOrder);
    assert(points == pointsSG, 'The number of points is invalid.');
  case 'tensor'
    [ nodes, weights, points ] = qd.constructTensorProduct(sdim, quadratureOrder);
    assert(points == pointsTP, 'The number of points is invalid.');
  otherwise
    error('The method type is unknown.');
  end

  grid = zeros(terms, points);

  norm = zeros(1, terms);

  for i = 1:terms
    f = Utils.toFunction(psi(i), x, 'rows');
    values = f(nodes);
    norm(i) = qd.computeNormalizationConstant(i, index);
    grid(i, :) = values .* weights / norm(i);
  end

  [ nodes, grid, norm ] = qd.finalize(sdim, nodes, grid, norm);

  assert(all(norm >= 0), 'Normalization constants cannot be negative.');
end
