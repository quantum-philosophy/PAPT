classdef GaussHermitePhysicists < Quadrature.Base
  methods
    function qd = GaussHermitePhysicists(x, psi, order, varargin)
      qd = qd@Quadrature.Base(x, psi, order, varargin{:});
    end
  end

  methods (Static)
    function [ nodes, weights ] = construct1D(order)
      %
      % `order' is the number of points involved.
      %
      [ nodes, weights ] = hermite_compute(order);

      nodes = nodes * sqrt(2);
      weights = weights / sqrt(pi);
    end

    function [ nodes, weights, points ] = constructTensorProduct(sdim, order)
      [ nodes1D, weights1D ] = Quadrature.GaussHermitePhysicists.construct1D(order);
      [ nodes, weights, points ] = ...
        Quadrature.constructTensorProduct(sdim, nodes1D, weights1D);
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order)
      %
      % order = 2^level - 1  --->  level = log2(order + 1)
      %
      level = ceil(log2(order + 1));

      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = nodes * sqrt(2);
      weights = weights / sqrt(pi^sdim);
    end
  end
end
