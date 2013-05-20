classdef Base < handle
  properties (SetAccess = 'private')
    expectation
    deviation
    transformation
    dimensionCount
  end

  methods
    function this = Base(varargin)
      options = Options(varargin{:});
      this.expectation = options.expectation;
      this.deviation = options.deviation;

      [ variance, correlation ] = this.correlate(options);
      variables = this.distribute(variance, correlation, options);

      this.transformation = ProbabilityTransformation.ReducedNormal( ...
        'variables', variables, 'threshold', options.threshold);

      this.dimensionCount = this.transformation.dimensionCount;
    end

    function data = evaluate(this, data)
      data = this.postprocess(this.transformation.evaluate(data));
    end

    function data = sample(this, sampleCount)
      data = this.postprocess(this.transformation.sample(sampleCount));
    end

    function display(this)
      options = Options( ...
        'Expectation', this.expectation, ...
        'Deviation', this.deviation, ...
        'Dimension', this.transformation.dimensionCount);
      display(options, 'Process variation');
    end
  end

  methods (Abstract, Access = 'protected')
    variables = distribute(this, variance, correlation, options)
  end

  methods (Access = 'private')
    [ variance, correlation ] = correlate(this, options);

    function data = postprocess(this, data)
      %
      % Join the local and global variations
      %
      data = bsxfun(@plus, data(:, end), data(:, 1:(end - 1)));

      data = this.expectation + this.deviation * data;
    end
  end
end
