function mc = constructMonteCarlo(method, c, samples)
  if nargin < 3, samples = c.samples; end

  %
  % Construct an appropriate solver.
  %
  hs = HotSpot.(method)(c.hotspotSet{:});
  f = @(rvs) hs.solve(c.dynamicPower, rvs);

  %
  % Construct a Monte Carlo sampler.
  %
  mc = MonteCarlo(f, [ hs.sdim, hs.cores, size(c.dynamicPower, 2) ], ...
    samples, c.stamp('prefix', method));
end
