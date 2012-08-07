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

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order, level)
      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = sqrt(2) * nodes;
      weights = weights / (sqrt(pi))^sdim;
    end

    function points = countSparseGridPoints(qd, sdim, order, level)
      points = sparse_grid_herm_size(sdim, level);
    end

    function norm = computeNormalizationConstant(qd, i, index)
      norm = prod(factorial(index(i, :) - 1));
    end
  end
end
