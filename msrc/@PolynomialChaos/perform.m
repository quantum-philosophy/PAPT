function [ E, C, out ] = perform(pc, f, samples)
  %
  % Output:
  %
  %   * E   - the expectation of `f',
  %   * C   - the covariance matrix of `f',
  %   * out - a bunch of samples.
  %
  %   NOTE: `out' is not used to compute the expectation and covariance.
  %

  if nargin < 3, samples = 10000; end

  sdim = pc.sdim;
  ddim = pc.ddim;

  %
  % Obtain the coefficients.
  %
  coeff = pc.gq.computeExpansion(f);

  %
  % Straight-forward stats.
  %
  E = coeff(:, 1);
  C = diag(sum(coeff(:, 2:end).^2 .* irep(pc.gq.norm(2:end), ddim, 1), 2));

  %
  % Now sampling.
  %
  out = transpose(pc.evaluateCustom(coeff, normrnd(0, 1, sdim, samples)));
end
