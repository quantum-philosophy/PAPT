function coeff = construct(pc, f, ddim)
  x = pc.x;
  psi = pc.psi;
  norm = pc.norm;
  terms = pc.terms;
  gq = pc.gq;

  coeff = zeros(ddim, terms);

  for i = 1:terms
    coeff(:, i) = gq.integrate(@(y) f(y) .* subs(psi(i), x, y), ddim) ./ norm(i);
  end
end
