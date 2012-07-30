classdef Probabilists < GaussianQuadrature.Base
  methods
    function gq = Probabilists(x, psi, order);
      gq = gq@GaussianQuadrature.Base(x, psi, order);
    end
  end

  methods (Static)
    [ nodes, weights, points ] = constructSparseGrid(sdim, order);
    [ nodes, weights, points ] = constructTensorProduct(sdim, order);
  end

  methods (Static, Access = 'protected')
    [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(x, psi, order);

    function [ nodes, plainGrid, niceGrid, norm ] = precomputeGrid(x, psi, order)
      %
      % A wrapper to cache the result of `doPrecomputeGrid'.
      %

      sdim = length(x);

      filename = [ 'PrQuadrature_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, plainGrid, niceGrid, norm ] = ...
          GaussianQuadrature.Probabilists.doPrecomputeGrid(x, psi, order);
        save(filename, 'nodes', 'plainGrid', 'niceGrid', 'norm');
      end
    end
  end
end
