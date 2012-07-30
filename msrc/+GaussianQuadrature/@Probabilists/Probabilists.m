classdef Probabilists < handle
  properties (SetAccess = 'protected')
    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid.
    %
    plainGrid

    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid multiplied by the
    % corresponding weight.
    %
    niceGrid

    %
    % The normalization coefficients of the expansion, i.e., <psi_i^2>.
    %
    norm
  end

  properties (SetAccess = 'protected')
    %
    % The evaluation points for the integration.
    %
    nodes

    %
    % The number of points.
    %
    points
  end

  methods
    function gq = Probabilists(x, psi, order);
      [ gq.nodes, gq.plainGrid, gq.niceGrid, gq.norm ] = gq.precomputeGrid(x, psi, order);
      gq.points = size(gq.nodes, 2);
    end
  end

  methods (Static, Access = 'protected')
    [ nodes, plainGrid, niceGrid, norm ] = doPrecomputeGrid(x, psi, order);
    [ nodes, weights, points ] = constructSparseGrid(sdim, level);
    [ nodes, weights, points ] = constructTensorProduct(sdim, level);

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
