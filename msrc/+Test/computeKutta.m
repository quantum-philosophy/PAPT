function [ exp, var, raw, time ] = computeKutta(c)
  if nargin < 1, c = Test.config(); end

  hs = HotSpot.Kutta(c.floorplan, c.hotspotConfig, c.hotspotLine);

  [ exp, var, raw, time ] = MonteCarlo.perform3D( ...
    @(rvs) hs.solve(c.dynamicPower, rvs), [ hs.sdim, c.cores, c.steps ], ...
    c.samples, c.stamp('Kutta'));

  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
