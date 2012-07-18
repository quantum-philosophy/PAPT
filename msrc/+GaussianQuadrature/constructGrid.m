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
