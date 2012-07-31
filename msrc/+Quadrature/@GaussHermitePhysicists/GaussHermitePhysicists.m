classdef GaussHermitePhysicists < Quadrature.Base
  methods
    function qd = GaussHermitePhysicists(varargin)
      qd = qd@Quadrature.Base(varargin{:});
    end
  end

  methods (Static)
    function level = orderToLevel(order)
      %
      % Description:
      %
      %   Here we are trying to estimate an appropriate level of the sparse grid
      %   such that the grid can integrate polynomial of the given order `order' exactly.
      %
      % Relationship:
      %
      %   order = 2^(level + 1) - 1
      %   level = log2(order + 1) - 1
      %
      % Source:
      %
      %   http://people.sc.fsu.edu/~jburkardt/m_src/sandia_sparse/sandia_sparse.html
      %
      level = ceil(log2(order + 1) - 1);
    end

    function points = countSparseGridPoints(sdim, order)
      level = Quadrature.GaussHermitePhysicists.orderToLevel(order);
      points = sparse_grid_herm_size(sdim, level);
    end

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
      level = Quadrature.GaussHermitePhysicists.orderToLevel(order);

      points = sparse_grid_herm_size(sdim, level);
      [ weights, nodes ] = sparse_grid_herm(sdim, level, points);

      nodes = nodes * sqrt(2);
      weights = weights / sqrt(pi^sdim);
    end
  end
end
