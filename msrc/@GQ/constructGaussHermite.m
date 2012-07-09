function [ nodes, weights, points ] = constructGaussHermite(dimension, level)
  %
  % Description:
  %
  %   Construct a sparse grid based on the Smolyak algorithm with
  %   the 1D Gauss-Hermite rule inside.
  %
  %   NOTE: The physicists' type of the Hermite polynomials is assumed,
  %   hence, the weight function is exp(-x^2).
  %

  points = sparse_grid_herm_size(dimension, level);
  [ weights, nodes ] = sparse_grid_herm(dimension, level, points);
end
