function out = evaluateCustom(pc, coeff, rvs)
  points = size(rvs, 2);

  ddim = pc.ddim;

  out = zeros(ddim, points);

  evaluator = pc.evaluator;

  for i = 1:ddim
    out(i, :) = evaluator(rvs, coeff(i, :));
  end
end
