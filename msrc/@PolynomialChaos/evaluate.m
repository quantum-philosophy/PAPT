function out = evaluate(pc, coeff, points)
  [ sdim, count ] = size(points);
  ddim = size(coeff, 1);

  vars = {};

  for i = 1:sdim
    vars{i} = points(i, :);
  end

  if ddim == 1
    e = matlabFunction((sum(coeff .* pc.psi)));
    out = e(vars{:});
  else
    e = matlabFunction((sum(coeff .* repmat(pc.psi, ddim, 1), 2)));
    out = reshape(e(vars{:}), count, ddim);
  end
end
