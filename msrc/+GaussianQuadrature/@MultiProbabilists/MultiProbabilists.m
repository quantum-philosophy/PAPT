classdef MultiProbabilists < GaussianQuadrature.Probabilists
  properties (SetAccess = 'protected')
    %
    % The normalization coefficients of the expansion, i.e., <psi_i^2>.
    %
    norm

    %
    % The number of terms in the PC expansion.
    %
    terms

    %
    % The number of stochastic dimensions.
    %
    ddim
  end

  methods
    function gq = MultiProbabilists(x, psi, index, ddim)
      gq = gq@GaussianQuadrature.Probabilists();
      [ gq.nodes, gq.plainGrid, grid, gq.norm ] = gq.precomputeGrid(x, psi, index);
      gq.points = size(gq.nodes, 2);
      gq.terms = length(gq.norm);

      %
      % Precompute the nice grid for the given number of
      % deterministic dimensions.
      %
      gq.ddim = ddim;
      gq.niceGrid = zeros(ddim, gq.points, gq.terms);
      for i = 1:gq.terms
        gq.niceGrid(:, :, i) = irep(grid(i, :), ddim, 1);
      end
    end

    function result = integrateWithChaos(gq, f, ddim, c)
      error('Deprecated.');
    end

    function result = integrateChaosProduct(gq, c1, c2)
      error('Deprecated.');
    end

    function coeff = computeExpansion(gq, f)
      grid = gq.niceGrid;
      terms = gq.terms;

      coeff = zeros(gq.ddim, terms);

      samples = f(gq.nodes);

      for i = 1:terms
        coeff(:, i) = sum(samples .* grid(:, :, i), 2);
      end
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
