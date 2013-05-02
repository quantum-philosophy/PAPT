classdef Base < Temperature.Analytical.Base
  properties (SetAccess = 'protected')
    process
    chaos
  end

  methods
    function this = Base(varargin)
      options = Options(varargin{:});

      this = this@Temperature.Analytical.Base(options);

      this.process = ProcessVariation( ...
        'expectation', LeakagePower.Base.Lnom, ...
        'deviation', 0.05 * LeakagePower.Base.Lnom, options);

      this.chaos = PolynomialChaos.Hermite( ...
        'inputCount', this.process.dimensionCount, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5), ...
        Options(varargin{:}));
    end

    function [ Texp, output ] = compute(this, Pdyn, varargin)
      chaos = this.chaos;

      coefficients = chaos.expand(@(rvs) ...
        transpose(this.solve(Pdyn, transpose(rvs), varargin{:})));

      [ processorCount, stepCount ] = size(Pdyn);

      Texp = reshape(coefficients(1, :), processorCount, stepCount);

      if nargout < 2, return; end

      outputCount = processorCount * stepCount;

      output.Tvar = reshape(sum(coefficients(2:end, :).^2 .* ...
        Utils.replicate(chaos.norm(2:end), 1, outputCount), 1), ...
        processorCount, stepCount);

      output.coefficients = reshape(coefficients, chaos.termCount, ...
        processorCount, stepCount);
    end

    function Tdata = sample(this, coefficients, sampleCount)
      rvs = normrnd(0, 1, sampleCount, this.process.dimensionCount);
      Tdata = this.evaluate(coefficients, rvs);
    end

    function Tdata = evaluate(this, coefficients, rvs)
      Tdata = this.chaos.evaluateSet(rvs, coefficients);
    end

    function display(this)
      display@Temperature.Analytical.Base(this);
      display(this.process);
      display(this.chaos);
    end
  end

  methods (Abstract)
    [ T, output ] = solve(this, Pdyn, rvs, varargin)
  end
end
