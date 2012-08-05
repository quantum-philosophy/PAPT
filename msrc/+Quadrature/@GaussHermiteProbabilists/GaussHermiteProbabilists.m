classdef GaussHermiteProbabilists < Quadrature.Base
  methods
    function qd = GaussHermiteProbabilists(varargin);
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Static)
    function points = countSparseGridPoints(sdim, order)
      [ ~, weights ] = nwspgr('gqn', sdim, order);
      points = length(weights);
    end

    function [ nodes, weights ] = construct1D(order)
      %
      % `order' is the number of points involved.
      %
      [ nodes, weights ] = nwspgr('gqn', 1, order);
    end

    function [ nodes, weights, points ] = constructTensorProduct(sdim, order)
      [ nodes1D, weights1D ] = Quadrature.GaussHermiteProbabilists.construct1D(order);
      [ nodes, weights, points ] = ...
        Quadrature.constructTensorProduct(sdim, nodes1D, weights1D);
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order);
      [ nodes, weights ] = nwspgr('gqn', sdim, order);
      nodes = transpose(nodes);
      weights = transpose(weights);
      points = size(weights, 2);
    end
  end
end
