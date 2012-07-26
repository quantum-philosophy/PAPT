function out = evaluate(pc, coeff, rvs)
  ddim = pc.ddim;

  out = zeros(ddim, pc.points);

  evaluator = pc.evaluator;

  for i = 1:ddim
    out(i, :) = evaluator(rvs, coeff(i, :));
  end
end
