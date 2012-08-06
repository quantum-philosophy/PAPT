classdef GaussJacobi < Quadrature.Base
  properties (Access = 'private')
    alpha
    beta
    a
    b
  end

  methods
    function qd = GaussJacobi(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Access = 'protected')
    function method = prepare(qd, method)
      if nargin < 2, method = struct(); end

      qd.alpha = method.jacobiAlpha;
      qd.beta = method.jacobiBeta;
      qd.a = method.jacobiA;
      qd.b = method.jacobiB;
    end

    function [ nodes, weights ] = construct1D(qd, order)
      [ nodes, weights ] = jacobi_compute(order, qd.alpha, qd.beta);
    end

    function [ nodes, weights, points ] = constructSparseGrid(qd, sdim, order)
      level = qd.orderToLevel(order);

      [ nodes, weights, points ] = Quadrature.constructSandiaGrid(...
        sdim, level, 9, qd.alpha, qd.beta);
    end

    function [ nodes, grid, norm ] = finalize(qd, sdim, nodes, grid, norm)
      nodes = ((qd.b - qd.a) / 2) * nodes + (qd.b + qd.a) / 2;
      norm = norm / (2^(qd.alpha + qd.beta + 1) * beta(qd.alpha + 1, qd.beta + 1))^sdim;
    end

    function norm = computeNormalizationConstant(qd, i, index)
      n = index(i, :) - 1;

      alpha = qd.alpha;
      beta = qd.beta;

      one = 2^(alpha + beta + 1) ./ (2 * n + alpha + beta + 1);
      two = gamma(n + alpha + 1) * gamma(n + beta + 1) ./ ...
        (gamma(n + 1) .* gamma(n + alpha + beta + 1));

      norm = prod(one .* two);
    end

    function points = countSparseGridPoints(qd, sdim, order)
      level = qd.orderToLevel(order);
      points = Quadrature.countSandiaPoints(sdim, level, 9, qd.alpha, qd.beta);
    end

    function level = orderToLevel(qd, order)
      level = ceil(log2(order + 1) - 1);
    end
  end
end
