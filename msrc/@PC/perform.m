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

  %
  % Obtain the coefficients.
  %
  coeffs = pc.construct(f);

  %
  % Straight-forward stats.
  %
  E = coeffs(1);
  C = sum(coeffs(2:end).^2 .* pc.norm(2:end));

  %
  % Now sampling.
  %
  vars = {};

  for i = 1:pc.dimension
    vars{i} = normrnd(0, 1, samples, 1);
  end

  e = matlabFunction(simplify(expand(sum(coeffs .* pc.psi))));
  out = e(vars{:});
end
