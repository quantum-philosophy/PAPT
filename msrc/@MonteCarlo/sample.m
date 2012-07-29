function out = sample(f, dims, samples, stamp)
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

    h = ibar('Monte Carlo: sample %d out of %d.', samples, done);

    m = tic;
    for i = 1:left
      out(done + i, :, :) = f(rvs(:, i));
      increase(h);
    end
    t = t + toc(m);

    name = [ 'MonteCarlo_', stamp, '_mcs', num2str(samples), '.mat' ];
    save(Utils.resolvePath(name, 'cache'), 'out', 't');
  end
end
