classdef Kutta < HotSpot.Base
  properties (Access = 'private')
    At
    Bt
    Ct
  end

  methods
    function hs = Kutta(floorplan, config, line)
      hs = hs@HotSpot.Base(floorplan, config, line);

      hs.At = - diag(1 ./ hs.C) * hs.G;
      hs.Bt = diag(1 ./ hs.C) * ...
        [ diag(ones(1, hs.cores)); zeros(hs.nodes - hs.cores, hs.cores) ];
      hs.Ct = diag(1 ./ hs.C) * hs.G * ones(hs.nodes, 1) * hs.Tamb;
    end

    function TT = solve(hs, Pdyn, rvs)
      [ cores, steps ] = size(Pdyn);
      assert(cores == hs.cores, 'The power profile is invalid.')

      %
      % General shortcuts.
      %
      At = hs.At;
      Bt = hs.Bt;
      Ct = hs.Ct;
      dt = hs.dt;
      Tamb = hs.Tamb;

      time = zeros(1, steps + 1);

      for i = 1:(steps + 1)
        time(i) = (i - 1) * dt;
      end

      T0 = ones(hs.nodes, 1) * Tamb;

      %
      % Initialize the leakage model.
      %
      leak = Leakage.Polynomial(Tamb, Pdyn, hs.map);

      [ ~, TT ] = ode45(...
        @(t, T) Ct + At * T + ...
          Bt * (Pdyn(:, min(round(t / dt) + 1, steps)) ...
            + leak.performAtGiven(T(1:cores, :), rvs)), ...
        time, T0);

      TT = transpose(TT(2:end, 1:cores));
    end
  end
end
