function [ nodes, weights, points ] = constructSparseGrid(sdim, order)
  %
  % order = 2^level - 1  --->  level = log2(order + 1)
  %
  level = ceil(log2(order + 1));

  points = sparse_grid_herm_size(sdim, level);

  [ weights, nodes ] = sparse_grid_herm(sdim, level, points);
end
