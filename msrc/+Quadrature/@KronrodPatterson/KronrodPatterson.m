classdef KronrodPatterson < Quadrature.Base
  methods
    function qd = KronrodPatterson(varargin);
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Static)
    function points = countTensorProductPoints(sdim, order)
      points = Inf;
    end

    function points = countSparseGridPoints(sdim, order)
      [ ~, weights ] = nwspgr('kpn', sdim, order);
      points = length(weights);
    end

    function [ nodes, weights ] = construct1D(order)
      error('It is meant to be sparse.');
    end

    function [ nodes, weights, points ] = constructTensorProduct(sdim, order)
      error('It is meant to be sparse.');
    end

    function [ nodes, weights, points ] = constructSparseGrid(sdim, order);
      [ nodes, weights ] = nwspgr('kpn', sdim, order);
      nodes = transpose(nodes);
      weights = transpose(weights);
      points = size(weights, 2);
    end
  end
end
