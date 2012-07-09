function [ nodes, weights, points ] = constructSparseGrid(dimension, level)
  filename = [ 'SG_d', num2str(dimension), '_l', num2str(level), '.mat' ];
  filename = Utils.resolvePath(filename);

  if exist(filename, 'file')
    load(filename);
  else
    [ nodes, weights, points ] = GQ.constructGaussHermite(dimension, level);
    save(filename, 'nodes', 'weights', 'points');
  end
end
