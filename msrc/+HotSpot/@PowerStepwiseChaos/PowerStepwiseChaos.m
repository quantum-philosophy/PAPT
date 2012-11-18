classdef PowerStepwiseChaos < HotSpot.StepwiseChaos
  methods
    function this = PowerStepwiseChaos(varargin)
      this = this@HotSpot.StepwiseChaos(varargin{:});
    end

    function [ Texp, Tvar, Tcoeff, Pexp, Pvar, Pcoeff ] = ...
      computeWithLeakage(this, Pdyn, leakage)

      [ processorCount, stepCount ] = size(Pdyn);
      assert(processorCount == this.processorCount);

      PdynT = Pdyn.';

      Tamb = this.ambientTemperature;

      chaos = this.chaos;

      ET = this.ET;
      DT = this.DT;
      B = this.B;

      Lnom = this.Lnom;
      Ldev = this.Ldev;
      LmapT = this.LmapT;

      %
      % Here we are going to store the stochastic power and temperature.
      %
      Pcoeff = zeros(chaos.termCount, processorCount, stepCount);
      Tcoeff = zeros(chaos.termCount, processorCount, stepCount);

      sample = @(rvs) leakage.evaluate(Lnom + rvs * LmapT * Ldev, Tamb);

      %
      % Perform the PC expansion and obtain the coefficients of
      % the current power.
      %
      Pcoeff(:, :, 1) = chaos.expand(sample);

      %
      % Add the dynamic power of the first step to the mean.
      %
      Pcoeff(1, :, 1) = Pcoeff(1, :, 1) + PdynT(1, :);

      sample = @(rvs, T) leakage.evaluate(Lnom + rvs * LmapT * Ldev, T);

      %
      % The first step is special because we do not have any expansion yet,
      % and the (projected) temperature is assumed to be zero.
      %
      Xcoeff = Pcoeff(:, :, 1) * DT;

      for i = 2:stepCount
        %
        % Calculate the previous temperature coefficients
        % (it is a real temperature now, i.e., in Kelvin).
        %
        Tcoeff(:, :, i - 1) = Xcoeff * B;
        Tcoeff(1, :, i - 1) = Tcoeff(1, :, i - 1) + Tamb;

        %
        % Perform the PC expansion.
        %
        Pcoeff(:, :, i) = chaos.expand(sample, Tcoeff(:, :, i - 1));

        %
        % Add the dynamic power to the mean.
        %
        Pcoeff(1, :, i) = Pcoeff(1, :, i) + PdynT(i, :);

        %
        % Compute new coefficients for each of the terms
        % of the PC expansion.
        %
        Xcoeff = Xcoeff * ET + Pcoeff(:, :, i) * DT;
      end

      %
      % Do not forget about the last coefficients.
      %
      Tcoeff(:, :, stepCount) = Xcoeff * B;
      Tcoeff(1, :, stepCount) = Tcoeff(1, :, stepCount) + Tamb;

      %
      % Compute the expectation of temperature.
      %
      Texp = squeeze(Tcoeff(1, :, :));

      if nargout < 2, return; end

      %
      % Compute the variance of temperature.
      %
      Tvar = zeros(stepCount, processorCount);
      norm = Utils.replicate(chaos.norm(2:end), 1, processorCount);
      for i = 1:stepCount
        Tvar(i, :) = sum(Tcoeff(2:end, :, i).^2 .* norm, 1);
      end
      Tvar = Tvar.';

      if nargout < 4, return; end

      %
      % Compute the expectation of power.
      %
      Pexp = squeeze(Pcoeff(1, :, :));

      if nargout < 5, return; end

      %
      % Compute the variance of power.
      %
      Pvar = zeros(stepCount, processorCount);
      for i = 1:stepCount
        Pvar(i, :) = sum(Pcoeff(2:end, :, i).^2 .* norm, 1);
      end
      Pvar = Pvar.';
    end
  end
end
