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

    function [ Texp, Tvar, coefficients ] = compute(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      chaos = this.chaos;

      coefficients = chaos.expand(@(L) this.solve(Pdyn, leakage, L));

      outputCount = processorCount * stepCount;

      Texp = reshape(coefficients(1, :), processorCount, stepCount);
      Tvar = reshape(sum(coefficients(2:end, :).^2 .* ...
        Utils.replicate(chaos.norm(2:end), 1, outputCount), 1), ...
        processorCount, stepCount);
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

  methods (Access = 'private')
    function T = solve(this, Pdyn, leakage, rvs)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      leak = leakage.evaluate;

      E = this.E;
      D = this.D;
      BT = this.BT;
      Tamb = this.ambientTemperature;
      process = this.process;

      sampleCount = size(rvs, 1);
      L = process.expectation + process.deviation * process.mapping * transpose(rvs);

      range = 1:processorCount;
      T = zeros(processorCount * stepCount, sampleCount);

      X = D * bsxfun(@plus, Pdyn(:, 1), leak(L, Tamb));
      T(range, :) = BT * X + Tamb;

      for i = 2:stepCount
        X = E * X + D * bsxfun(@plus, Pdyn(:, i), leak(L, T(range, :)));
        range = range + processorCount;
        T(range, :) = BT * X + Tamb;
      end

      T = transpose(T);
    end
  end
end
