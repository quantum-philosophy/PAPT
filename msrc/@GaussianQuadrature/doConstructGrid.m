function [ nodes, weights, count ] = doConstructGrid(sdim, level)
  %
  % Description:
  %
  %   Construct a sparse grid based on the Smolyak algorithm with
  %   the 1D Gauss-Hermite rule inside.
  %
  %   NOTE: The physicists' type of the Hermite polynomials is assumed,
  %   hence, the weight function is exp(-x^2).
  %

  count = sparse_grid_herm_size(sdim, level);
  [ weights, nodes ] = sparse_grid_herm(sdim, level, count);
end
