classdef KronrodPatterson < Quadrature.Base
  methods
    function qd = KronrodPatterson(varargin);
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function [ nodes, weights ] = construct1D(qd, order)
      error('It is meant to be sparse.');
    end

    function [ nodes, weights, points ] = constructTensorProduct(qd, sdim, order)
      error('It is meant to be sparse.');
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order, level);
      [ nodes, weights ] = nwspgr('kpn', sdim, order);
      nodes = transpose(nodes);
      weights = transpose(weights);
      points = size(weights, 2);
    end
  end

  methods (Static)
    function points = countTensorProductPoints(sdim, order)
      points = Inf;
    end

    function points = countSparseGridPoints(sdim, order, level)
      [ ~, weights ] = nwspgr('kpn', sdim, order);
      points = length(weights);
    end
  end
end
