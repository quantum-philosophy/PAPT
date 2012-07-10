function coeff = construct(pc, f, axes)
  x = pc.x;
  psi = pc.psi;
  norm = pc.norm;
  count = pc.count;
  gq = pc.gq;

  coeff = zeros(axes, count);

  for i = 1:count
    coeff(:, i) = gq.integrate(@(y) f(y) .* subs(psi(i), x, y), axes) ./ norm(i);
  end
end
