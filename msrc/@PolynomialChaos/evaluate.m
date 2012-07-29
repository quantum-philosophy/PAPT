function result = evaluate(pc, coeff, rvs)
  [ ddim, terms ] = size(coeff);
  assert(ddim == pc.ddim, 'The deterministic dimension is invalid.');
  assert(terms == pc.terms, 'The number of terms is invalid.');

  mterms = size(pc.rvProd, 1);
  points = size(rvs, 2);

  rvPower = pc.rvPower;
  rvProd = zeros(mterms, points);

  for i = 1:mterms
    rvProd(i, :) = prod(realpow(rvs, irep(rvPower(:, i), 1, points)), 1);
  end

  result = (coeff * pc.coeffMap) * rvProd;
end
