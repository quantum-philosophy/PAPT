function [ E, V ] = perform3D(f, dims, samples)
  %
  % Description:
  %
  %   Samples the given function `f' a number of times
  %   equal to `samples', and each time it is expected that
  %   the function will return its execution trace (in time).
  %   For the each moment of time, the expectation and variance
  %   are computed.
  %
  % Output:
  %
  %   * E   - the profile of expectations of `f',
  %   * V   - the profile of variances of `f'.
  %

  if nargin < 2, dims = [ 1, 1, 1 ]; end
  if nargin < 3, samples = 10000; end

  %
  % Stochastic dimension.
  %
  if numel(dims) < 1
    sdim = 1;
  else
    sdim = dims(1);
  end

  %
  % Deterministic dimension.
  %
  if numel(dims) < 2
    ddim = 1;
  else
    ddim = dims(2);
  end

  %
  % Time dimension.
  %
  if numel(dims) < 3
    tdim = 1;
  else
    tdim = dims(3);
  end

  %
  % First, we sample.
  %
  rvs = normrnd(0, 1, sdim, samples);
  out = zeros(samples, ddim, tdim);

  h = waitbar(0, sprintf('Monte Carlo sampling: %d/%d.', 0, samples));

  for i = 1:samples
    out(i, :, :) = f(rvs(:, i));
    waitbar(i / samples, h, sprintf('Monte Carlo sampling: %d/%d.', i, samples));
  end

  close(h);

  %
  % Compute the expectation and variance.
  %
  E = zeros(ddim, tdim);
  V = zeros(ddim, tdim);

  for i = 1:tdim
    snapshot = out(:, :, i);
    E(:, i) = mean(snapshot);
    V(:, i) = diag(cov(snapshot));
  end
end
