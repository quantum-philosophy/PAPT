classdef MonteCarlo < handle
  properties (Access = 'private')
    data
  end

  methods
    function mc = MonteCarlo(varargin)
      mc.data = MonteCarlo.sample(varargin{:});
    end

    function [ exp, var, raw, I ] = sampleSequential(mc, count, I)
      %
      % Description:
      %
      %   Takes the loaded samples always from the very beginning.
      %   `I' specifies the samples to exclude.
      %

      [ totalCount, ddim, tdim ] = size(mc.data);

      assert(count <= totalCount, 'The number of samples is too large.');

      if nargin < 3, I = []; end

      % All the samples.
      J = 1:totalCount;
      % Remove the excluded samples.
      J(I) = [];
      % Take the needed amount from the rest.
      J = J(1:count);

      raw = mc.data(J, :, :);

      % Accumulate the banned samples.
      I = [ I J ];

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

    function [ exp, var, raw, I ] = sampleAccumulative(mc, count, I)
      %
      % Description:
      %
      %   Takes the loaded samples at random, but never repeats.
      %   `I' maintains already drawn samples.
      %

      [ totalCount, ddim, tdim ] = size(mc.data);

      assert(count <= totalCount, 'The number of samples is too large.');

      if nargin < 3, I = []; end

      % Determine how many more to take.
      more = max(0, count - length(I));

      if more > 0
        % All the samples.
        J = 1:totalCount;
        % Remove already taken.
        J(I) = [];
        % Merge.
        I = [ I J(randperm(length(J), more)) ];
      end

      raw = mc.data(I(1:count), :, :);

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

    function [ exp, var, raw, I ] = sampleWithReplacement(mc, count, I)
      %
      % Description:
      %
      %   Takes the loaded samples at random, possibly with repetitions.
      %   `I' specifies the range of samples to choose from.
      %

      [ totalCount, ddim, tdim ] = size(mc.data);

      if nargin < 3
        I = randi(totalCount, count, 1);
      end

      raw = mc.data(I, :, :);
      raw = raw(randi(count, count, 1), :, :);

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
    out = sample(f, dims, samples, stamp)
  end
end
