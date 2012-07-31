classdef Physicists < GaussianQuadrature.Base
  methods
    function gq = Physicists(x, psi, order, varargin)
      gq = gq@GaussianQuadrature.Base(x, psi, order, varargin{:});
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
      [ nodes1D, weights1D ] = GaussianQuadrature.Physicists.construct1D(order);
      [ nodes, weights, points ] = ...
        GaussianQuadrature.constructTensorProduct(sdim, nodes1D, weights1D);
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order)
      %
      % order = 2^level - 1  --->  level = log2(order + 1)
      %
      level = ceil(log2(order + 1));

      error('Not clear how to choose the level.');

      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = nodes * sqrt(2);
      weights = weights / sqrt(pi^sdim);
    end
  end
end
