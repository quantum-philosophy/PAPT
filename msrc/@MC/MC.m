classdef MC < handle
  methods (Static)
    function [ mu, cm, out ] = perform(f, dimension, samples)
      if nargin < 2, dimension = 1; end
      if nargin < 3, samples = 10000; end

      rvs = normrnd(0, 1, samples, dimension);
      out = zeros(samples, 1);

      for i = 1:samples
        out(i, :) = f(rvs(i, :));
      end

      mu = mean(out);
      cm = cov(out);
    end
  end
end
