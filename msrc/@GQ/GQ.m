classdef GQ < handle
  %
  % (*) The quadrature rule that is used here (Gauss-Hermite) assumes
  % the physicists' weight function. We need the probabilists's one,
  % therefore, a proper conversion is to be performed.
  %
  properties (SetAccess = 'private')
    %
    % The dimension of the grid, i.e., the number of variables.
    %
    dimension

    %
    % The level of the grid.
    %
    level

    %
    % The total number of nodes.
    %
    points

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
    function gq = GQ(dimension, level)
      if nargin < 2, level = 2; end

      gq.dimension = dimension;
      gq.level = level;

      [ nodes, gq.weights, gq.points ] = GQ.constructSparseGrid(dimension, level);

      %
      % See (*) to justify the need of sqrt(2).
      %
      gq.nodes = nodes * sqrt(2);
    end

    function result = integrate(gq, f, axes)
      if nargin < 3, axes = 1; end

      points = gq.points;
      nodes = gq.nodes;

      samples = zeros(axes, points);

      for i = 1:points
        samples(:, i) = f(nodes(:, i));
      end

      %
      % See (*) to understand why we do not need 2 next to pi here.
      %
      result = sum(samples .* repmat(gq.weights, axes, 1), 2) ./ ...
        pi^(gq.dimension / 2);
    end
  end

  methods (Static, Access = 'private')
    [ nodes, weights, points ] = constructSparseGrid(dimension, level);
    [ nodes, weights, points ] = constructGaussHermite(dimension, level);
  end
end
