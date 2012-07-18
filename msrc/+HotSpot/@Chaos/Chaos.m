classdef Chaos < HotSpot.Analytic
  properties (Access = 'private')
    %
    % The polynomial chaos.
    %
    pc
  end

  methods
    function hs = Chaos(floorplan, config, line, order)
      hs = hs@HotSpot.Analytic(floorplan, config, line);

      if nargin < 4, order = 2; end

      %
      % Initialize the PC expansion.
      %
      hs.pc = PolynomialChaos(hs.sdim, order);
    end

    function [ ExpT, VarT ] = solve(hs, Pdyn)
      [ cores, steps ] = size(Pdyn);
      assert(cores == hs.cores, 'The power profile is invalid.')

      %
      % General shortcuts.
      %
      nodes = hs.nodes;
      E = hs.E;
      D = hs.D;
      BT = hs.BT;
      Tamb = hs.Tamb;

      %
      % Shortcuts for the PC expansion.
      %
      pc = hs.pc;
      terms = pc.terms;

      %
      % Here we are going to store the stochastic temperature.
      %
      trace = ChaosTrace(cores, terms, steps);

      %
      % Initialize the leakage model.
      %
      leak = Leakage.Polynomial(Tamb, Pdyn, hs.map, trace, pc);

      %
      % Perform the PC expansion and obtain the coefficients of
      % the current power.
      %
      Pcoeff = pc.construct(@leak.performAtAmbient, cores);

      %
      % Add the dynamic power of the first step to the mean.
      %
      Pcoeff(:, 1) = Pcoeff(:, 1) + Pdyn(:, 1);

      %
      % The storage of the current (projected) temperature coefficients.
      %
      Tcoeff = zeros(nodes, terms);

      %
      % The first step is special because we do not have any expansion yet,
      % and the (projected) temperature is assumed to be zero.
      %
      for i = 1:terms
        Tcoeff(:, i) = D * Pcoeff(:, i);
      end

      for i = 2:steps
        %
        % Calculate the previous temperature coefficients
        % (it is real temperature now, i.e., in Kelvin).
        %
        trace.coeff(:, :, i - 1) = BT * Tcoeff;
        trace.coeff(:, 1, i - 1) = trace.coeff(:, 1, i - 1) + Tamb;

        %
        % Perform the PC expansion.
        %
        leak.advance();
        Pcoeff = pc.construct(@leak.performAtCurrent, cores);

        %
        % Add the dynamic power to the mean.
        %
        Pcoeff(:, 1) = Pcoeff(:, 1) + Pdyn(:, i);

        %
        % Compute new coefficients for each of the terms
        % of the PC expansion.
        %
        for j = 1:terms
          Tcoeff(:, j) = E * Tcoeff(:, j) + D * Pcoeff(:, j);
        end
      end

      %
      % Do not forget about the last coefficients.
      %
      trace.coeff(:, :, steps) = BT * Tcoeff;
      trace.coeff(:, 1, steps) = trace.coeff(:, 1, steps) + Tamb;

      %
      % Compute the expectation.
      %
      ExpT = squeeze(trace.coeff(:, 1, :));

      %
      % Compute the variance.
      %
      % NOTE: The correlation matrix is full of zeros.
      %
      VarT = zeros(cores, steps);
      norm = irep(pc.norm(2:end), cores, 1);
      for i = 1:steps
        VarT(:, i) = sum(trace.coeff(:, 2:end, i).^2 .* norm, 2);
      end
    end
  end
end
