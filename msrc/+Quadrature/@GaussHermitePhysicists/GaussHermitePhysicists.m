classdef GaussHermitePhysicists < Quadrature.Base
  methods
    function qd = GaussHermitePhysicists(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function [ nodes, weights ] = construct1D(qd, order)
      [ nodes, weights ] = hermite_compute(order);

      nodes = sqrt(2) * nodes;
      weights = weights / sqrt(pi);
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order)
      level = qd.orderToLevel(order);

      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = sqrt(2) * nodes;
      weights = weights / (sqrt(pi))^sdim;
    end

    function [ nodes, grid, norm ] = finalize(qd, sdim, nodes, grid, norm)
    end

    function norm = computeNormalizationConstant(qd, i, index)
      norm = prod(factorial(index(i, :) - 1));
    end

    function level = orderToLevel(qd, order)
      %
      % Relationship:
      %
      %   order = 2^(level + 1) - 1
      %   level = log2(order + 1) - 1
      %
      level = ceil(log2(order + 1) - 1);
    end

    function points = countSparseGridPoints(qd, sdim, order)
      level = qd.orderToLevel(order);
      points = sparse_grid_herm_size(sdim, level);
    end
  end
end
