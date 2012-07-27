function [ E, V, out, t ] = perform3D(f, dims, samples, stamp)
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
  if nargin < 4, stamp = []; end

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

  out = zeros(0, ddim, tdim);
  done = 0;
  t = 0;

  if ~isempty(stamp)
    pattern = [ 'MonteCarlo_', stamp, '_mcs(.+).mat' ];
    match = Utils.findCache(pattern);

    if ~isempty(match)
      load(match);
      done = min(samples, size(out, 1));
    end
  end

  left = samples - done;
  out = [ out(1:done, :, :); zeros(left, ddim, tdim) ];

  if left > 0
    rvs = normrnd(0, 1, sdim, left);

    h = waitbar(0, sprintf('Monte Carlo sampling: %d/%d.', done, samples));

    m = tic;
    for i = 1:left
      out(done + i, :, :) = f(rvs(:, i));
      waitbar((done + i) / samples, h, ...
        sprintf('Monte Carlo sampling: %d/%d.', done + i, samples));
    end
    t = t + toc(m);

    close(h);

    name = [ 'MonteCarlo_', stamp, '_mcs', num2str(samples), '.mat' ];
    save(Utils.resolvePath(name, 'cache'), 'out', 't');
  end

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
