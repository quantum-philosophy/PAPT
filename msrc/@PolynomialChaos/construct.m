function coeff = construct(pc, f, ddim)
  cq = pc.cq;
  terms = pc.terms;

  coeff = zeros(ddim, terms);

  for i = 1:terms
    coeff(:, i) = cq.integrateWithChaos(f, ddim, i);
  end

  coeff = coeff ./ irep(pc.norm, ddim, 1);
end
