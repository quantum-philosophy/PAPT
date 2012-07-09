function coeffs = construct(pc, f)
  x = pc.x;
  psi = pc.psi;
  norm = pc.norm;
  count = pc.count;
  gq = pc.gq;

  coeffs = zeros(1, count);

  for i = 1:count
    coeffs(i) = gq.integrate(@(y) f(y) * subs(psi(i), x, y)) / norm(i);
  end
end
