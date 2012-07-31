classdef KronrodPatterson < Quadrature.Base
  methods
    function qd = KronrodPatterson(x, psi, order, varargin);
      qd = qd@Quadrature.Base(x, psi, order, varargin{:});
    end
  end

  methods (Access = 'protected')
    function points = countTensorProductPoints(qd, sdim, order)
      points = Inf;
    end
  end

  methods (Static)
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
