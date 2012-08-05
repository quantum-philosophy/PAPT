function [ nodes, weights, points ] = constructTensorProduct(qd, sdim, order)
  [ nodes1D, weights1D ] = qd.construct1D(order);

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
  assert(points == order^sdim);
end
