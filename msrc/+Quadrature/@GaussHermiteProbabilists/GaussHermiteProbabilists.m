classdef GaussHermiteProbabilists < Quadrature.Base
  methods
    function qd = GaussHermiteProbabilists(varargin);
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function [ nodes, weights ] = construct1D(qd, order)
      %
      % `order' is the number of points involved.
      %
      [ nodes, weights ] = nwspgr('gqn', 1, order);
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order);
      [ nodes, weights ] = nwspgr('gqn', sdim, order);
      nodes = transpose(nodes);
      weights = transpose(weights);
      points = size(weights, 2);
    end

    function norm = computeNormalizationConstant(qd, i, index)
      norm = prod(factorial(index(i, :) - 1));
    end

    function points = countSparseGridPoints(qd, sdim, order)
      [ ~, weights ] = nwspgr('gqn', sdim, order);
      points = length(weights);
    end
  end
end
