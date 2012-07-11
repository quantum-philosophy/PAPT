function out = evaluate(pc, coeff, points)
  [ sdim, count ] = size(points);
  ddim = size(coeff, 1);

  out = zeros(ddim, count);

  for i = 1:ddim
    out(i, :) = pc.evaluator(points, coeff(i, :));
  end
end
