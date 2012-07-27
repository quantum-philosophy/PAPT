function result = evaluateCustom(pc, coeff, rvs)
  mterms = size(pc.rvProd, 1);
  points = size(rvs, 2);

  rvPower = pc.rvPower;
  rvProd = zeros(mterms, points);

  for i = 1:mterms
    rvProd(i, :) = prod(realpow(rvs, irep(rvPower(:, i), 1, points)), 1);
  end

  result = (coeff * pc.coeffMap) * rvProd;
end
