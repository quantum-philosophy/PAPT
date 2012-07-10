function [ E, C, out ] = perform(pc, f, dims, count)
  %
  % Output:
  %
  %   * E   - the expectation of `f',
  %   * C   - the covariance matrix of `f',
  %   * out - a bunch of samples.
  %
  %   NOTE: `out' is not used to compute the expectation and covariance.
  %

  if nargin < 4, count = 10000; end

  sdim = dims(1);
  ddim = dims(2);

  if sdim ~= pc.dimension
    error('The dimensions do not match each other.');
  end

  %
  % Obtain the coefficients.
  %
  coeff = pc.construct(f, ddim);

  %
  % Straight-forward stats.
  %
  E = coeff(:, 1);
  C = diag(sum(coeff(:, 2:end).^2 .* repmat(pc.norm(2:end), ddim, 1), 2));

  %
  % Now sampling.
  %
  out = transpose(pc.evaluate(coeff, normrnd(0, 1, sdim, count)));
end
