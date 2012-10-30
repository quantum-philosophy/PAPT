classdef Chaos < HotSpot.Analytic
  properties (SetAccess = 'private')
    options

    Lnom
    rvMap
    rvCount
  end

  methods
    function this = Chaos(floorplan, config, line, varargin)
      this = this@HotSpot.Analytic(floorplan, config, line);

      this.Lnom = LeakagePower.Lnom;

      [ this.rvMap, this.rvCount ] = ...
        ProcessVariation.analyze(floorplan, this.Lnom);

      this.options = Options( ...
        'inputCount', this.rvCount, ...
        'outputCount', 0, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5));
    end

    function [ Texp, Tvar ] = computeWithLeakage(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      chaos = PolynomialChaos.Hermite(@(L) this.solve(Pdyn, leakage, L), ...
        this.options, 'outputCount', processorCount * stepCount);

      Texp = reshape(chaos.expectation, processorCount, stepCount);
      Tvar = reshape(chaos.variance, processorCount, stepCount);
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
