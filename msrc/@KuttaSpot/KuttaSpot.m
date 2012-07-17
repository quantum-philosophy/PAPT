classdef KuttaSpot < HotSpot
  properties (Access = 'private')
    At
    Bt
    Ct
  end

  methods
    function ks = KuttaSpot(floorplan, hsConfig, hsLine)
      ks = ks@HotSpot(floorplan, hsConfig, hsLine);

      ks.At = - diag(1 ./ ks.C) * ks.G;
      ks.Bt = diag(1 ./ ks.C) * ...
        [ diag(ones(1, ks.cores)); zeros(ks.nodes - ks.cores, ks.cores) ];
      ks.Ct = diag(1 ./ ks.C) * ks.G * ones(ks.nodes, 1) * ks.Tamb;
    end

    function TT = solve(ks, Pdyn, rvs)
      [ cores, steps ] = size(Pdyn);
      assert(cores == ks.cores, 'The power profile is invalid.')

      %
      % General shortcuts.
      %
      At = ks.At;
      Bt = ks.Bt;
      Ct = ks.Ct;
      dt = ks.dt;
      Tamb = ks.Tamb;

      time = zeros(1, steps + 1);

      for i = 1:(steps + 1)
        time(i) = (i - 1) * dt;
      end

      T0 = ones(ks.nodes, 1) * Tamb;

      %
      % Initialize the leakage model.
      %
      leak = Leakage(Tamb, Pdyn, ks.map);

      [ ~, TT ] = ode45(...
        @(t, T) Ct + At * T + ...
          Bt * (Pdyn(:, min(round(t / dt) + 1, steps)) ...
            + leak.performAtGiven(T(1:cores, :), rvs)), ...
        time, T0);

      TT = transpose(TT(2:end, 1:cores));
    end
  end
end
