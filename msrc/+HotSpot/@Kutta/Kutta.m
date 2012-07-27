classdef Kutta < HotSpot.Base
  properties (Access = 'private')
    At
    Bt
  end

  methods
    function hs = Kutta(floorplan, config, line)
      hs = hs@HotSpot.Base(floorplan, config, line);

      hs.At = - diag(1 ./ hs.C) * hs.G;
      hs.Bt = diag(1 ./ hs.C) * ...
        [ diag(ones(1, hs.cores)); zeros(hs.nodes - hs.cores, hs.cores) ];
    end

    function TT = solve(hs, Pdyn, rvs)
      [ cores, steps ] = size(Pdyn);
      assert(cores == hs.cores, 'The power profile is invalid.')

      %
      % General shortcuts.
      %
      At = hs.At;
      Bt = hs.Bt;
      dt = hs.dt;
      Tamb = hs.Tamb;

      %
      % Initialize the leakage model.
      %
      leak = Leakage.Polynomial(Tamb, Pdyn, hs.pca);

      TT = zeros(steps, cores);

      T0 = ones(1, hs.nodes) * Tamb;

      for i = 1:steps
        [ ~, T0 ] = ode45(...
          @(t, T) At * (T - Tamb) + Bt * (Pdyn(:, i) ...
              + leak.performAtGiven(rvs, T(1:cores))), ...
          [ 0, dt ], T0);

          T0 = T0(end, :);
          TT(i, :) = T0(1:cores);
      end

      TT = transpose(TT);
    end
  end
end
