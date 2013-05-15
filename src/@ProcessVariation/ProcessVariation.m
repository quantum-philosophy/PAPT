classdef ProcessVariation < handle
  properties (SetAccess = 'protected')
    expectation
    deviation
    mapping
    dimensionCount
  end

  methods
    function this = ProcessVariation(varargin)
      options = Options(varargin{:});

      this.expectation = options.expectation;
      this.deviation = options.deviation;

      this.construct(options);
    end

    function [ u, n, z ] = sample(this, count)
      z = randn(this.dimensionCount, count);
      n = this.mapping * z;
      u = this.expectation + this.deviation * n;
    end

    function display(this)
      options = Options( ...
        'Expectation', this.expectation, ...
        'Deviation', this.deviation, ...
        'Dimension', this.dimensionCount);
      display(options, 'Process variation');
    end
  end

  methods (Access = 'private')
    mapping = construct(this, options)
  end
end
