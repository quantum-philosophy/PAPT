classdef GaussHermitePhysicists < Quadrature.Base
  properties Constant
    nodeScale = sqrt(2);
    weightScale = 1 / sqrt(pi);
  end

  methods
    function qd = GaussHermitePhysicists(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function [ nodes, weights ] = construct1D(qd, order)
      %
      % `order' is the number of points involved.
      %
      [ nodes, weights ] = hermite_compute(order);

      nodes = qd.nodeScale * nodes;
      weights = qd.weightScale * weights;
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order)
      level = Quadrature.GaussHermitePhysicists.orderToLevel(order);

      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = qd.nodeScale * nodes;
      weights = qd.weightScale^sdim * weights;
    end

    function norm = computeNorm(qd, i, index)
      norm = prod(factorial(index(i, :) - 1));
    end
  end

  methods (Static)
    function level = orderToLevel(order)
      %
      % Relationship:
      %
      %   order = 2^(level + 1) - 1
      %   level = log2(order + 1) - 1
      %
      level = ceil(log2(order + 1) - 1);
    end

    function points = countSparseGridPoints(sdim, order)
      level = Quadrature.GaussHermitePhysicists.orderToLevel(order);
      points = sparse_grid_herm_size(sdim, level);
    end
  end
end
