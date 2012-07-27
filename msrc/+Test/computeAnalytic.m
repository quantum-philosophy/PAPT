function [ exp, var, raw, time ] = computeAnalytic(c)
  if nargin < 1, c = Test.config(); end

  hs = HotSpot.Analytic(c.floorplan, c.hotspotConfig, c.hotspotLine);

  [ exp, var, raw, time ] = MonteCarlo.perform3D( ...
    @(rvs) hs.solve(c.dynamicPower, rvs), [ hs.sdim, c.cores, c.steps ], ...
    c.samples, c.stamp('Analytic'));

  exp = Utils.toCelsius(exp);
  raw = Utils.toCelsius(raw);
end
