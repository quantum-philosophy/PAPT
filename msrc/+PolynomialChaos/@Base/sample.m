function [ exp, var, raw ] = sample(pc, f, points)
  %
  % Output:
  %
  %   * exp - the expectation of `f',
  %   * var - the variance of `f',
  %   * raw - a bunch of samples.
  %
  %   NOTE: `out' is not used to compute the expectation and covariance.
  %

  if nargin < 3, points = 10000; end

  sdim = pc.sdim;
  ddim = pc.ddim;

  %
  % Obtain the coefficients.
  %
  coeff = pc.computeExpansion(f);

  %
  % Straight-forward stats.
  %
  exp = coeff(:, 1);
  var = diag(sum(coeff(:, 2:end).^2 .* irep(pc.norm(2:end), ddim, 1), 2));

  %
  % Now sampling.
  %
  raw = transpose(pc.evaluate(coeff, pc.generateSampleNodes(points)));
end
