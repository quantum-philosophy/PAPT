classdef ChaosQuadrature < handle
  properties (Access = 'private')
    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid.
    %
    plainGrid

    %
    % Precomputed value of each of the polynomials in the PC expansion
    % in each of the points of the sparse grid multiplied by the
    % corresponding weight and divided by the corresponding normalization
    % coefficient.
    %
    niceGrid
  end

  properties (SetAccess = 'private')
    %
    % The evaluation points for the integration.
    %
    nodes

    %
    % The number of points.
    %
    count
  end

  methods
    function cq = ChaosQuadrature(x, psi, order)
      [ cq.nodes, cq.plainGrid, cq.niceGrid ] = ...
        cq.precomputeGrid(x, psi, order);
      cq.count = size(cq.nodes, 2);
    end

    function result = integrateWithChaos(cq, f, ddim, c)
      samples = f(cq.nodes);
      result = sum(samples .* irep(cq.niceGrid(c, :), ddim, 1), 2);
    end

    function result = integrateChaosProduct(cq, c1, c2)
      result = sum(cq.plainGrid(c1, :) .* cq.niceGrid(c2, :));
    end
  end

  methods (Static, Access = 'private')
    [ nodes, plainGrid, niceGrid ] = doPrecomputeGrid(x, psi, order);

    function [ nodes, plainGrid, niceGrid ] = precomputeGrid(x, psi, order)
      %
      % A wrapper to cache the result of `doPrecomputeGrid'.
      %

      sdim = length(x);

      filename = [ 'CQ_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, plainGrid, niceGrid ] = ...
          ChaosQuadrature.doPrecomputeGrid(x, psi, order);
        save(filename, 'nodes', 'plainGrid', 'niceGrid');
      end
    end
  end
end
