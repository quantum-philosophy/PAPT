function [ nodes, weights, points ] = constructSparseGrid(sdim, level);
  [ nodes, weights ] = nwspgr('gqn', sdim, level);
  nodes = transpose(nodes);
  weights = transpose(weights);

  points = size(weights, 2);
end
