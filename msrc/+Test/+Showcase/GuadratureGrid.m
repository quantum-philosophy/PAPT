init;

sdim = 2;
order = 10;

requiredOrder = order + 1;

fprintf('Random variables:    %d\n', sdim);
fprintf('Polynomial order:    %d\n', order);
fprintf('Integration order:   %d\n', requiredOrder);

[ nodes, weights, points ] = ...
  GaussianQuadrature.Physicists.constructSparseGrid(sdim, requiredOrder);
fprintf('Physicists sparse:   %d (%d negative)\n', points, length(find(weights < 0)));

[ nodes, weights, points ] = ...
  GaussianQuadrature.Probabilists.constructSparseGrid(sdim, requiredOrder);
fprintf('Probabilists sparse: %d (%d negative)\n', points, length(find(weights < 0)));

[ nodes, weights, points ] = ...
  GaussianQuadrature.Physicists.constructTensorProduct(sdim, requiredOrder);
fprintf('Physicists tensor:   %d (%d negative)\n', points, length(find(weights < 0)));

[ nodes, weights, points ] = ...
  GaussianQuadrature.Probabilists.constructTensorProduct(sdim, requiredOrder);
fprintf('Probabilists tensor: %d (%d negative)\n', points, length(find(weights < 0)));
