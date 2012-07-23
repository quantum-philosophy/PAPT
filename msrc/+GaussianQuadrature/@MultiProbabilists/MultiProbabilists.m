classdef MultiProbabilists < GaussianQuadrature.Probabilists
  properties (SetAccess = 'protected')
    %
    % The normalization coefficients of the expansion, i.e., <psi_i^2>.
    %
    norm
  end

  methods
    function gq = MultiProbabilists(x, psi, index)
      gq = gq@GaussianQuadrature.Probabilists();
      [ gq.nodes, gq.plainGrid, gq.niceGrid, gq.norm ] = gq.precomputeGrid(x, psi, index);
      gq.points = size(gq.nodes, 2);
    end

    function result = integrateWithChaos(gq, f, ddim, c)
      error('Deprecated.');
    end

    function result = integrateChaosProduct(gq, c1, c2)
      error('Deprecated.');
    end

    function result = integrateWithNormalizedChaos(gq, f, ddim, c)
      samples = f(gq.nodes);
      result = sum(samples .* irep(gq.niceGrid(c, :), ddim, 1), 2);
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
