classdef Chaos < HotSpot.Analytic & ProcessVariation
  properties (Access = 'protected')
    chaos
  end

  methods
    function this = Chaos(floorplan, config, line, varargin)
      this = this@HotSpot.Analytic(floorplan, config, line);
      this = this@ProcessVariation(floorplan, varargin{:});

      this.chaos = PolynomialChaos.Hermite( ...
        'inputCount', this.rvCount, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5), ...
        Options(varargin{:}));
    end

    function [ Texp, Tvar, coefficients ] = ...
      computeWithLeakage(this, Pdyn, leakage)

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
      samples = normrnd(0, 1, sampleCount, this.rvCount);
      Tdata = this.chaos.evaluateSet(samples, coefficients);
    end
  end

  methods (Access = 'private')
    function T = solve(this, Pdyn, leakage, L)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      E = this.E;
      D = this.D;
      BT = this.BT;
      Tamb = this.ambientTemperature;

      sampleCount = size(L, 1);
      L = this.Lnom + this.rvMap * transpose(L);

      range = 1:processorCount;
      T = zeros(processorCount * stepCount, sampleCount);

      X = D * bsxfun(@plus, Pdyn(:, 1), leakage.evaluate(L, Tamb));
      T(range, :) = BT * X + Tamb;

      for i = 2:stepCount
        X = E * X + D * bsxfun(@plus, Pdyn(:, i), leakage.evaluate(L, T(range, :)));
        range = range + processorCount;
        T(range, :) = BT * X + Tamb;
      end

      T = transpose(T);
    end
  end
end
