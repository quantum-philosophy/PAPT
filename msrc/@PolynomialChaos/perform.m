function [ E, C, out ] = perform(pc, f, dims, samples)
  %
  % Output:
  %
  %   * E   - the expectation of `f',
  %   * C   - the covariance matrix of `f',
  %   * out - a bunch of samples.
  %
  %   NOTE: `out' is not used to compute the expectation and covariance.
  %

  if nargin < 4, samples = 10000; end

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

  vars = {};

  for i = 1:sdim
    vars{i} = normrnd(0, 1, samples, 1);
  end

  %
  % Now sampling.
  %
  if ddim == 1
    %
    % Just slightly simplified construction for
    % the single-space-dimension case.
    %
    e = matlabFunction(sum(coeff .* pc.psi));
    out = e(vars{:});
  else
    e = matlabFunction(sum(coeff .* repmat(pc.psi, ddim, 1), 2));
    out = reshape(e(vars{:}), samples, ddim);
  end
end
