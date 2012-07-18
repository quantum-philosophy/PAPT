function [ nodes, weights, points ] = constructSparseGrid(sdim, level)
  points = sparse_grid_herm_size(sdim, level);
  [ weights, nodes ] = sparse_grid_herm(sdim, level, points);
end
