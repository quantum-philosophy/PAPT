function mc = sample(method, c, samples)
  if nargin < 2, c = Test.config(); end
  if nargin < 3, samples = c.samples; end

  stamp = c.stamp('prefix', method, 'samples', samples);

  hs = HotSpot.(method)(c.hotspotSet{:});
  f = @(rvs) hs.solve(c.dynamicPower, rvs);

  mc = MonteCarlo(f, [ hs.sdim, c.cores, c.steps ], samples, stamp);
end
