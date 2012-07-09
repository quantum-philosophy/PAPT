classdef MC < handle
  methods (Static)
    function [ E, C, out ] = perform(f, dimension, samples)
      %
      % Output:
      %
      %   * E   - the expectation of `f',
      %   * C   - the covariance matrix of `f',
      %   * out - the samples used to compute the stats.
      %

      if nargin < 2, dimension = 1; end
      if nargin < 3, samples = 10000; end

      rvs = normrnd(0, 1, samples, dimension);
      out = zeros(samples, 1);

      for i = 1:samples
        out(i, :) = f(rvs(i, :));
      end

      E = mean(out);
      C = cov(out);
    end
  end
end
