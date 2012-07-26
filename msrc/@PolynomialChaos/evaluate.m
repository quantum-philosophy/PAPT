function out = evaluate(pc, coeff, points)
  [ sdim, count ] = size(points);
  ddim = size(coeff, 1);
  evaluator = pc.evaluator;

  out = zeros(ddim, count);

  for i = 1:ddim
    out(i, :) = evaluator(points, coeff(i, :));
  end
end
