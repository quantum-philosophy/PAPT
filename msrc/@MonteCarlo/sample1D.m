function [ E, C, out ] = sample1D(f, dims, samples)
  %
  % Description:
  %
  %   Samples the given function `f' a number of times
  %   equal to `samples' and outputs its expectation,
  %   covariance, and the computed samples.
  %
  % Output:
  %
  %   * E   - the expectation of `f',
  %   * C   - the covariance matrix of `f',
  %   * out - the samples used to compute the stats.
  %

  if nargin < 3, samples = 10000; end

  sdim = dims(1);
  ddim = dims(2);

  rvs = normrnd(0, 1, sdim, samples);
  out = zeros(ddim, samples);

  for i = 1:samples
    out(:, i) = f(rvs(:, i));
  end

  out = transpose(out);

  E = mean(out);
  C = cov(out);
end
