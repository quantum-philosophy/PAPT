function [ nodes, weights, points ] = constructSparseGrid(sdim, order);
  [ nodes, weights ] = nwspgr('gqn', sdim, order);
  nodes = transpose(nodes);
  weights = transpose(weights);

  points = size(weights, 2);
end
