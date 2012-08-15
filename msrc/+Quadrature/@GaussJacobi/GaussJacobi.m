classdef GaussJacobi < Quadrature.Base
  properties (Access = 'private')
    alpha
    beta
    a
    b
  end

  methods
    function qd = GaussJacobi(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function method = prepare(qd, method)
      method = prepare@Quadrature.Base(qd, method);

      qd.alpha = method.jacobiAlpha;
      qd.beta = method.jacobiBeta;
      qd.a = method.jacobiA;
      qd.b = method.jacobiB;
    end

    function [ nodes, weights ] = construct1D(qd, order)
      [ nodes, weights ] = jacobi_compute(order, qd.alpha, qd.beta);
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order, level)
      o = ones(1, sdim);
      rule  = o * 9;
      alpha = o * qd.alpha;
      beta  = o * qd.beta;
      tol = sqrt(eps);

      total_points = sparse_grid_mixed_size_total(sdim, level, rule);

      points = sparse_grid_mixed_size(sdim, level, rule, alpha, beta, tol);

      sparse_unique_index = sparse_grid_mixed_unique_index( ...
        sdim, level, rule, alpha, beta, tol, points, total_points);

      [ sparse_order, sparse_index ] = sparse_grid_mixed_index( ...
        sdim, level, rule, points, total_points, sparse_unique_index);

      nodes = sparse_grid_mixed_point(sdim, level, rule, ...
        alpha, beta, points, sparse_order, sparse_index);

      weights = sparse_grid_mixed_weight(sdim, level, rule, ...
        alpha, beta, points, total_points, sparse_unique_index);
    end

    function [ nodes, grid, norm ] = finalize(qd, sdim, nodes, grid, norm)
      nodes = ((qd.b - qd.a) / 2) * nodes + (qd.b + qd.a) / 2;
      norm = norm / (2^(qd.alpha + qd.beta + 1) * beta(qd.alpha + 1, qd.beta + 1))^sdim;
    end

    function norm = computeNormalizationConstant(qd, i, index)
      n = index(i, :) - 1;

      alpha = qd.alpha;
      beta = qd.beta;

      one = 2^(alpha + beta + 1) ./ (2 * n + alpha + beta + 1);
      two = gamma(n + alpha + 1) .* gamma(n + beta + 1) ./ ...
        (gamma(n + 1) .* gamma(n + alpha + beta + 1));

      norm = prod(one .* two);
    end
  end

  methods (Static)
    function points = countSparseGridPoints(sdim, order, level)
      o = ones(1, sdim);
      rule  = o * 9;
      alpha = o * 69;
      beta  = o * 69;
      tol = sqrt(eps);

      points = sparse_grid_mixed_size(sdim, level, rule, alpha, beta, tol);
    end
  end
end
