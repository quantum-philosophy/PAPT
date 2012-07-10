function [ E, C ] = perform3D(f, dims, samples)
  if nargin < 3, samples = 10000; end

  sdim = dims(1);
  ddim = dims(2);
  tdim = dims(3);

  rvs = normrnd(0, 1, sdim, samples);
  out = zeros(samples, ddim, tdim);

  for i = 1:samples
    out(i, :, :) = f(rvs(:, i));
  end

  E = zeros(ddim, tdim);
  C = zeros(ddim, tdim);

  for i = 1:tdim
    shot = out(:, :, i);
    E(:, i) = mean(shot);
    C(:, i) = diag(cov(shot));
  end
end
