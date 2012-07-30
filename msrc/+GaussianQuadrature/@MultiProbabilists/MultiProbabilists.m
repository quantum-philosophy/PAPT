classdef MultiProbabilists < GaussianQuadrature.Probabilists
  methods
    function gq = MultiProbabilists(x, psi, index)
      gq = gq@GaussianQuadrature.Probabilists(x, psi, index);
    end
  end

  methods (Static, Access = 'protected')
    [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(x, psi, index);

    function [ nodes, plainGrid, niceGrid, norm ] = precomputeGrid(x, psi, index)
      %
      % A wrapper to cache the result of `doPrecomputeGrid'.
      %

      sdim = length(x);
      order = PolynomialChaos.indexToOrder(index);

      filename = [ 'MultiPrQuadrature_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, plainGrid, niceGrid, norm ] = ...
          GaussianQuadrature.MultiProbabilists.doPrecomputeGrid(x, psi, index);
        save(filename, 'nodes', 'plainGrid', 'niceGrid', 'norm');
      end
    end
  end
end
