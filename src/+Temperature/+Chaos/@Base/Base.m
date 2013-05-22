classdef Base < Temperature.Analytical.Base
  properties (SetAccess = 'protected')
    process
    chaos
  end

  methods
    function this = Base(varargin)
      options = Options(varargin{:});

      this = this@Temperature.Analytical.Base(options.temperatureOptions);
      this.process = ProcessVariation.(options.processModel)( ...
        options.processOptions);
      this.chaos = PolynomialChaos.Hermite( ...
        'inputCount', this.process.dimensionCount, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5), ...
        options.chaosOptions);
    end

    function [ Texp, output ] = compute(this, Pdyn, varargin)
      chaos = this.chaos;
      process = this.process;

      coefficients = chaos.expand(@(rvs) transpose(this.solve( ...
        Pdyn, transpose(process.evaluate(rvs)), varargin{:})));

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
      Tdata = this.chaos.sample(sampleCount, coefficients);
    end

    function Tdata = evaluate(this, coefficients, rvs)
      Tdata = this.chaos.evaluate(rvs, coefficients);
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
