classdef Chaos < HotSpot.Analytic
  properties (SetAccess = 'protected')
    process
    chaos
  end

  methods
    function this = Chaos(varargin)
      options = Options(varargin{:});

      this = this@HotSpot.Analytic(options);

      this.process = ProcessVariation( ...
        'expectation', LeakagePower.Lnom, ...
        'deviation', 0.05 * LeakagePower.Lnom, options);

      this.chaos = PolynomialChaos.Hermite( ...
        'inputCount', this.process.dimensionCount, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5), ...
        Options(varargin{:}));
    end

    function [ Texp, Tvar, coefficients ] = compute(this, Pdyn, varargin)
      chaos = this.chaos;

      coefficients = chaos.expand(@(rvs) this.solve(Pdyn, rvs, varargin{:}));

      [ processorCount, stepCount ] = size(Pdyn);

      Texp = reshape(coefficients(1, :), processorCount, stepCount);

      if nargout < 2, return; end

      outputCount = processorCount * stepCount;

      Tvar = reshape(sum(coefficients(2:end, :).^2 .* ...
        Utils.replicate(chaos.norm(2:end), 1, outputCount), 1), ...
        processorCount, stepCount);

      if nargout < 3, return; end

      coefficients = reshape(coefficients, chaos.termCount, ...
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
      display@HotSpot.Analytic(this);
      display(this.process);
      display(this.chaos);
    end
  end

  methods (Access = 'protected')
    function [ T, P ] = solve(this, Pdyn, rvs, varargin)
      options = Options(varargin{:});

      process = this.process;
      options.L = process.expectation + ...
        process.deviation * process.mapping * transpose(rvs);

      [ T, P ] = this.([ 'solve', options.method ])(Pdyn, options);

      T = transpose(T);
      P = transpose(P);
    end
  end
end
