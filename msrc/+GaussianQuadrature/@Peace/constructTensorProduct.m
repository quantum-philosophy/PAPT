function [ nodes, weights, points ] = constructTensorProduct(sdim, level);
  [ nodes1D, weights1D ] = nwspgr('gqn', 1, level);

  nodes = {};
  weights = {};

  for i = 1:sdim
    nodes{end + 1} = nodes1D;
    weights{end + 1} = weights1D;
  end

  [ nodes, weights ] = tensor_product(nodes, weights);
  nodes = transpose(nodes);
  weights = transpose(weights);

  points = size(weights, 2);
  assert(points == level^sdim);
end
