function [ nodes, weights, points ] = constructSandiaGrid(sdim, level, rule, alpha, beta)
  one = ones(sdim, 1);

  rule = one * rule;
  alpha = one * alpha;
  beta = one * beta;
  tol = sqrt(eps);

  point_total_num = sparse_grid_mixed_size_total(sdim, level, rule);

  points = sparse_grid_mixed_size(sdim, level, rule, alpha, beta, tol);

  sparse_unique_index = sparse_grid_mixed_unique_index( ...
    sdim, level, rule, alpha, beta, tol, points, point_total_num);

  [ sparse_order, sparse_index ] = sparse_grid_mixed_index( ...
    sdim, level, rule, points, point_total_num, sparse_unique_index);

  nodes = sparse_grid_mixed_point(sdim, level, rule, ...
    alpha, beta, points, sparse_order, sparse_index);

  weights = sparse_grid_mixed_weight(sdim, level, rule, ...
    alpha, beta, points, point_total_num, sparse_unique_index);
end
