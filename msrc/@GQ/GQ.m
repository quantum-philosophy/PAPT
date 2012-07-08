classdef GQ < handle
  methods (Static)
    function result = integrate(f, dimension, level)
      if nargin < 3, level = 2; end

      [ nodes, weights, points ] = GQ.getSparseGrid(dimension, level);

      nodes = nodes * sqrt(2);

      samples = zeros(1, points);
      for i = 1:points
        samples(i) = f(nodes(:, i));
      end

      result = sum(samples .* weights) / pi^(dimension / 2);
    end
  end

  methods (Static, Access = 'private')
    function [ nodes, weights, points ] = getSparseGrid(dimension, level)
      filename = [ 'SG_d', num2str(dimension), '_l', num2str(level), '.mat' ];
      filename = Utils.resolvePath(filename);

      if exist(filename, 'file')
        load(filename);
      else
        [ nodes, weights, points ] = GQ.constructGaussHermite(dimension, level);
        save(filename, 'nodes', 'weights', 'points');
      end
    end

    function [ nodes, weights, points ] = constructGaussHermite(dimension, level)
      points = sparse_grid_herm_size(dimension, level);
      [ weights, nodes ] = sparse_grid_herm(dimension, level, points);
    end
  end
end
