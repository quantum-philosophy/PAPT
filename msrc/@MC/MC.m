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

      out = f(normrnd(0, 1, samples, dimension));

      E = mean(out);
      C = cov(out);
    end
  end
end
