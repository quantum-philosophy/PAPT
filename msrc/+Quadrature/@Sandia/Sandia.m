classdef Sandia < Quadrature.Base
  methods
    function qd = Sandia(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Static)
    function points = countTensorProductPoints(sdim, order)
      points = Inf;
    end

    function points = countSparseGridPoints(sdim, order)
      rule = ones(sdim, 1) * 5;
      level = Quadrature.Sandia.orderToLevel(order);
      points = sparse_grid_mixed_size(sdim, level, rule, 0, 0, sqrt(eps));
    end

    function level = orderToLevel(order)
      level = ceil(log2(order + 1) - 1);
    end

    function [ nodes, weights ] = construct1D(order)
      error('Is not meant to be used.');
    end

    function [ nodes, weights, points ] = constructTensorProduct(sdim, order)
      error('Is not meant to be used.');
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order)
      alpha = 0;
      beta = 0;
      rule = ones(sdim, 1) * 5;
      tol = sqrt(eps);

      level = Quadrature.Sandia.orderToLevel(order);

      point_total_num = sparse_grid_mixed_size_total ( sdim, level, rule );

      points = sparse_grid_mixed_size ( sdim, level, rule, alpha, beta, tol );

      sparse_unique_index = sparse_grid_mixed_unique_index ( ...
        sdim, level, rule, alpha, beta, tol, points, point_total_num );

      [ sparse_order, sparse_index ] = sparse_grid_mixed_index ( ...
        sdim, level, rule, points, point_total_num, sparse_unique_index );

      nodes = sparse_grid_mixed_point ( sdim, level, rule, ...
        alpha, beta, points, sparse_order, sparse_index );

      weights = sparse_grid_mixed_weight ( sdim, level, rule, ...
        alpha, beta, points, point_total_num, sparse_unique_index );

      nodes = nodes * sqrt(2);
      weights = weights / sqrt(pi^sdim);
    end
  end
end
