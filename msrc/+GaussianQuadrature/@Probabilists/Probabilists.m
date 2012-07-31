classdef Probabilists < GaussianQuadrature.Base
  methods
    function gq = Probabilists(x, psi, order, varargin);
      gq = gq@GaussianQuadrature.Base(x, psi, order, varargin{:});
    end
  end

  methods (Static)
    function [ nodes, weights ] = construct1D(order)
      %
      % `order' is the number of points involved.
      %
      [ nodes, weights ] = nwspgr('gqn', 1, order);
    end

    function [ nodes, weights, points ] = constructTensorProduct(sdim, order)
      [ nodes1D, weights1D ] = GaussianQuadrature.Probabilists.construct1D(order);
      [ nodes, weights, points ] = ...
        GaussianQuadrature.constructTensorProduct(sdim, nodes1D, weights1D);
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order);
      [ nodes, weights ] = nwspgr('gqn', sdim, order);
      nodes = transpose(nodes);
      weights = transpose(weights);
      points = size(weights, 2);
    end
  end
end
