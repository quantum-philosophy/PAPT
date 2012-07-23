function coeff = construct(pc, f, ddim)
  qd = pc.qd;
  terms = pc.terms;

  coeff = zeros(ddim, terms);

  for i = 1:terms
    coeff(:, i) = qd.integrateWithNormalizedChaos(f, ddim, i);
  end
end
