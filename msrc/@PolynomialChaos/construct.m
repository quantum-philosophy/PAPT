function coeff = construct(pc, f, ddim)
  qd = pc.qd;
  terms = pc.terms;

  coeff = zeros(ddim, terms);

  for i = 1:terms
    coeff(:, i) = qd.integrateWithChaos(f, ddim, i);
  end

  coeff = coeff ./ irep(pc.norm, ddim, 1);
end
