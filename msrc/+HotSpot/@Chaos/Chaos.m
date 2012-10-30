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

      inputCount = this.rvCount;
      outputCount = processorCount * stepCount;

      E = this.E;
      D = this.D;
      BT = this.BT;
      Tamb = this.ambientTemperature;
      Lnom = this.Lnom;
      rvMap = this.rvMap;

      zeros = @uninit;

      function data = solve(samples)
        sampleCount = size(samples, 1);
        samples = Lnom + rvMap * transpose(samples);

        data = zeros(sampleCount, outputCount);

        for i = 1:sampleCount
          L = samples(:, i);
          T = zeros(processorCount, stepCount);

          X = D * (Pdyn(:, 1) + leakage.evaluate(L, Tamb));
          T(:, 1) = BT * X + Tamb;

          for j = 2:stepCount
            X = E * X + D * (Pdyn(:, j) + leakage.evaluate(L, T(:, j - 1)));
            T(:, j) = BT * X + Tamb;
          end

          data(i, :) = reshape(T, 1, []);
        end
      end

      options = this.options;
      options.outputCount = outputCount;

      chaos = PolynomialChaos.Hermite(@solve, this.options);

      display(chaos);

      Texp = reshape(chaos.expectation, processorCount, stepCount);
      Tvar = reshape(chaos.variance, processorCount, stepCount);
    end
  end
end
