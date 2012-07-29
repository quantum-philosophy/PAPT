classdef MonteCarlo < handle
  properties (Access = 'private')
    data
    taken
  end

  methods
    function mc = MonteCarlo(f, dims, count, stamp)
      if nargin < 2, dims = [ 1, 1, 1 ]; end
      if nargin < 3, count = 10000; end
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
          done = min(count, size(out, 1));
        end
      end

      left = count - done;
      out = [ out(1:done, :, :); zeros(left, ddim, tdim) ];

      if left > 0
        rvs = normrnd(0, 1, sdim, left);

        h = ibar('Monte Carlo: sample %d out of %d.', count, done);

        m = tic;
        for i = 1:left
          out(done + i, :, :) = f(rvs(:, i));
          increase(h);
        end
        t = t + toc(m);

        name = [ 'MonteCarlo_', stamp, '_mcs', num2str(count), '.mat' ];
        save(Utils.resolvePath(name, 'cache'), 'out', 't');
      end

      mc.taken = [];
      mc.data = out;
    end

    function [ exp, var, raw ] = sample(mc, count)
      [ totalCount, ddim, tdim ] = size(mc.data);

      assert(count <= totalCount, 'The number of samples is too large.');

      takenCount = length(mc.taken);

      if takenCount < count
        extraCount = count - takenCount;

        ids = 1:totalCount;
        ids(mc.taken) = [];

        leftCount = length(ids);

        extraTaken = randperm(leftCount, extraCount);

        mc.taken = [ mc.taken extraTaken ];
      end

      raw = mc.data(mc.taken(1:count), :, :);

      %
      % Compute the expectation and variance.
      %
      exp = zeros(ddim, tdim);
      var = zeros(ddim, tdim);

      for i = 1:tdim
        snapshot = raw(:, :, i);
        exp(:, i) = mean(snapshot);
        var(:, i) = diag(cov(snapshot));
      end
    end
  end

  methods (Static)
    [ E, C, out ] = sample1D(f, dims, samples);
  end
end
