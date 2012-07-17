classdef GaussianQuadrature < handle
  %
  % (*) The quadrature rule that is used here (Gauss-Hermite) assumes
  % the physicists' weight function. We need the probabilists' one,
  % therefore, a proper conversion is to be performed.
  %
  properties (SetAccess = 'private')
    %
    % The dimension of the grid, i.e., the number of variables.
    %
    sdim

    %
    % The level of the grid.
    %
    level

    %
    % The total number of nodes.
    %
    count

    %
    % The nodes.
    %
    nodes

    %
    % The corresponding weights.
    %
    weights
  end

  methods
    function gq = GaussianQuadrature(sdim, order)
      gq.sdim = sdim;
      gq.order = order;

      [ nodes, gq.weights, gq.count ] = ...
        GaussianQuadrature.constructGrid(sdim, order);

      %
      % See (*) to justify the need of sqrt(2).
      %
      gq.nodes = nodes * sqrt(2);
    end

    function result = integrate(gq, f, ddim)
      count = gq.count;
      nodes = gq.nodes;

      samples = zeros(ddim, count);

      for i = 1:count
        samples(:, i) = f(nodes(:, i));
      end

      %
      % See (*) to understand why we do not need 2 next to pi here.
      %
      result = sum(samples .* irep(gq.weights, ddim, 1), 2) ./ pi^(gq.sdim / 2);
    end
  end

  methods (Static)
    function [ nodes, weights, count ] = constructGrid(sdim, order)
      %
      % A wrapper to cache the result of `doConstructGrid'.
      %

      filename = [ 'QG_d', num2str(sdim), '_o', num2str(order), '.mat' ];
      filename = Utils.resolvePath(filename, 'cache');

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, weights, count ] = ...
          GaussianQuadrature.doConstructGrid(sdim, order);
        save(filename, 'nodes', 'weights', 'count');
      end
    end

    function level = computeLevel(order)
      %
      % Description:
      %
      %   An n-point Gaussian quadrature rule is exact for
      %   polynomials of order (2 * n - 1) or less.
      %
      % Inputs:
      %
      %   * order - the order of the polynomial.
      %

      level = ceil((order + 1) / 2);
    end
  end

  methods (Static, Access = 'private')
    [ nodes, weights, count ] = doConstructGrid(sdim, level);
  end
end
