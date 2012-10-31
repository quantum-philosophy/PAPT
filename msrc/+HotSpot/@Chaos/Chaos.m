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
      chaos.expandPermanent(@(L) this.solve(Pdyn, leakage, L));

      Texp = reshape(chaos.expectation, processorCount, stepCount);
      Tvar = reshape(chaos.variance, processorCount, stepCount);
      coefficients = reshape(chaos.coefficients, chaos.termCount, ...
        processorCount, stepCount);
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
