classdef ChaosHeat < HotSpot
  properties (Constant)
    %
    % The portion of the leakage power in terms of the dynamic one.
    %
    leakageRatioInDynamic = 1.0;

    %
    % The temperature for the peak leakage above.
    %
    leakageTemperature = Utils.toKelvin(100);
  end

  properties (Access = 'private')
    %
    % The stochastic dimension of the analysis, i.e.,
    % the number of uncertain parameters.
    %
    params

    %
    % The order of the polynomial chaos (PC) expansion.
    %
    order

    %
    % The `params'-dimensional polynomial chaos.
    %
    pc
  end

  methods
    function ch = ChaosHeat(floorplan, hsConfig, hsLine, order)
      %
      % Inputs:
      %
      %   * floorplan, hsConfig, hsLine - the same as for HotSpot,
      %
      %   * order - the order of the PC expansion.
      %

      ch = ch@HotSpot(floorplan, hsConfig, hsLine);

      if nargin < 4, order = 2; end

      %
      % Initialize the PC expansion. One core corresponds to
      % one random variable for now.
      %
      ch.params = ch.cores;
      ch.order = order;
      ch.pc = PolynomialChaos(ch.params, ch.order);
    end

    function [ ExpT, VarT ] = solve(ch, Pdyn)
      [ cores, steps ] = size(Pdyn);
      if cores ~= ch.cores, error('The power profile is invalid.'); end

      %
      % General shortcuts.
      %
      nodes = ch.nodes;
      E = ch.E;
      D = ch.D;

      %
      % Shortcuts for the PC expansion.
      %
      pc = ch.pc;
      terms = ch.pc.terms;

      %
      % Initialize the leakage power model.
      %
      % NOTE: It is computed for each core separately.
      %
      leakage = Leakage(ChaosHeat.leakageRatioInDynamic * mean(Pdyn, 2), ...
        ChaosHeat.leakageTemperature);

      %
      % Here we are going to store the stochastic temperature.
      %
      Tcoeff = zeros(nodes, terms, steps);

      %
      % Perform the PC expansion.
      %
      Pcoeff = zeros(cores, terms);

      %
      % Add the dynamic power of the first step to the mean.
      %
      Pcoeff(:, 1) = Pcoeff(:, 1) + Pdyn(:, 1);

      %
      % The first step is special because we do not have any expansion yet,
      % and the (projected) temperature is assumed to be zero.
      %
      for i = 1:terms
        Tcoeff(:, i, 1) = D * Pcoeff(:, i);
      end

      for i = 2:steps
        %
        % Perform the PC expansion.
        %
        Pcoeff = zeros(cores, terms);

        %
        % Add the dynamic power to the mean.
        %
        Pcoeff(:, 1) = Pcoeff(:, 1) + Pdyn(:, i);

        %
        % Compute new coefficients for each of the terms
        % of the PC expansion.
        %
        for j = 1:terms
          Tcoeff(:, j, i) = E * Tcoeff(:, j, i - 1) + D * Pcoeff(:, j);
        end
      end

      %
      % Compute the expectation.
      %
      ExpT = ch.BT * squeeze(Tcoeff(:, 1, :)) + ch.Tamb;

      %
      % Compute the variance.
      %
      % NOTE: We do not need the covariance matrix here since
      % there are no correlations between cores after the PC expansion.
      %
      VarT = zeros(cores, steps);
      norm = repmat(pc.norm(2:end), cores, 1);
      for i = 1:steps
        VarT(:, i) = sum((ch.BT * Tcoeff(:, 2:end, i)).^2 .* norm, 2);
      end
    end
  end
end
