classdef StepwiseChaos < HotSpot.Chaos
  properties (Access = 'private')
    ET
    DT
    B
    rvMapT
  end

  methods
    function this = StepwiseChaos(varargin)
      this = this@HotSpot.Chaos(varargin{:});

      this.ET = this.E.';
      this.DT = this.D.';
      this.B = this.BT.';
      this.rvMapT = this.rvMap.';
    end

    function [ Texp, Tvar, coefficients ] = ...
      computeWithLeakage(this, Pdyn, leakage)

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

    function Tdata = sample(this, coefficients, sampleCount)
      samples = normrnd(0, 1, sampleCount, this.rvCount);
      Tdata = this.chaos.evaluateSet(samples, coefficients);
    end
  end
end
