classdef Chaos < HotSpot.Analytic
  properties (Access = 'private')
    Lnom

    rvMap
    rvCount

    chaos

    %
    % For the stepwise version.
    %
    ET
    DT
    B
    rvMapT
  end

  methods
    function this = Chaos(floorplan, config, line, varargin)
      this = this@HotSpot.Analytic(floorplan, config, line);

      this.Lnom = LeakagePower.Lnom;

      [ this.rvMap, this.rvCount ] = ...
        ProcessVariation.analyze(floorplan, this.Lnom);

      this.chaos = PolynomialChaos.Hermite( ...
        'inputCount', this.rvCount, ...
        'order', 4, ...
        'quadratureOptions', Options( ...
          'method', 'tensor', ...
          'order', 5));

      this.ET = this.E.';
      this.DT = this.D.';
      this.B = this.BT.';
      this.rvMapT = this.rvMap.';
    end

    function [ Texp, Tvar, coefficients ] = computeWithLeakage(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      chaos = this.chaos;
      chaos.expandPermanent(@(L) this.solve(Pdyn, leakage, L));

      Texp = reshape(chaos.expectation, processorCount, stepCount);
      Tvar = reshape(chaos.variance, processorCount, stepCount);
      coefficients = reshape(chaos.coefficients, chaos.termCount, processorCount, stepCount);
    end

    function [ Texp, Tvar, coefficients ] = computeWithLeakageStepwise(this, Pdyn, leakage)
      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      PdynT = Pdyn.';

      Tamb = this.ambientTemperature;

      Lnom = this.Lnom;

      chaos = this.chaos;

      ET = this.ET;
      DT = this.DT;
      B = this.B;
      rvMapT = this.rvMapT;

      %
      % Here we are going to store the stochastic temperature.
      %
      coefficients = zeros(chaos.termCount, processorCount, stepCount);

      sample = @(L) leakage.evaluate(Lnom + L * rvMapT, Tamb);

      %
      % Perform the PC expansion and obtain the coefficients of
      % the current power.
      %
      Pcoeff = chaos.expand(sample);

      %
      % Add the dynamic power of the first step to the mean.
      %
      Pcoeff(1, :) = Pcoeff(1, :) + PdynT(1, :);

      sample = @(L, T) leakage.evaluate(Lnom + L * rvMapT, T);

      %
      % The first step is special because we do not have any expansion yet,
      % and the (projected) temperature is assumed to be zero.
      %
      Tcoeff = Pcoeff * DT;

      for i = 2:stepCount
        %
        % Calculate the previous temperature coefficients
        % (it is a real temperature now, i.e., in Kelvin).
        %
        coefficients(:, :, i - 1) = Tcoeff * B;
        coefficients(1, :, i - 1) = coefficients(1, :, i - 1) + Tamb;

        %
        % Perform the PC expansion.
        %
        Pcoeff = chaos.expand(sample, coefficients(:, :, i - 1));

        %
        % Add the dynamic power to the mean.
        %
        Pcoeff(1, :) = Pcoeff(1, :) + PdynT(i, :);

        %
        % Compute new coefficients for each of the terms
        % of the PC expansion.
        %
        Tcoeff = Tcoeff * ET + Pcoeff * DT;
      end

      %
      % Do not forget about the last coefficients.
      %
      coefficients(:, :, stepCount) = Tcoeff * B;
      coefficients(1, :, stepCount) = coefficients(1, :, stepCount) + Tamb;

      %
      % Compute the expectation.
      %
      Texp = squeeze(coefficients(1, :, :));

      %
      % Compute the variance.
      %
      Tvar = zeros(stepCount, processorCount);
      norm = Utils.replicate(chaos.norm(2:end), 1, processorCount);
      for i = 1:stepCount
        Tvar(i, :) = sum(coefficients(2:end, :, i).^2 .* norm, 1);
      end
      Tvar = Tvar.';
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
