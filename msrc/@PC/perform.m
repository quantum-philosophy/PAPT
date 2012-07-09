function [ mu, sigma2 ] = perform(pc, f)
  coeffs = pc.construct(f);

  mu = coeffs(1);
  sigma2 = sum(coeffs.^2 .* pc.norm);
end
