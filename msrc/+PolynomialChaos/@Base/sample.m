function [ E, C, out ] = sample(pc, f, points)
  %
  % Output:
  %
  %   * E   - the expectation of `f',
  %   * C   - the covariance matrix of `f',
  %   * out - a bunch of samples.
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
  E = coeff(:, 1);
  C = diag(sum(coeff(:, 2:end).^2 .* irep(pc.norm(2:end), ddim, 1), 2));

  %
  % Now sampling.
  %
  out = transpose(pc.evaluate(coeff, pc.generateSampleNodes(points)));
end
