function points = countSandiaPoints(sdim, level, rule, alpha, beta)
  one = ones(sdim, 1);

  rule = one * rule;
  alpha = one * alpha;
  beta = one * beta;
  tol = sqrt(eps);

  points = sparse_grid_mixed_size(sdim, level, rule, alpha, beta, tol);
end
