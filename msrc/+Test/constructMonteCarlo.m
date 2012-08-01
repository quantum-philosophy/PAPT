function mc = constructMonteCarlo(c)
  %
  % Construct an appropriate solver.
  %
  hs = HotSpot.(c.assessmentMethod)(c.hotspotArguments{:}, 'none');

  assert(hs.sdim == (PrincipalComponent.globalCount + hs.cores), ...
    'The number of stochastic dimensions is invalid.');

  f = @(rvs) hs.solve(c.dynamicPower, rvs);

  %
  % Construct a Monte Carlo sampler.
  %
  mc = MonteCarlo(f, [ hs.sdim, c.cores, c.steps ], ...
    c.monteCarloTotalSamples, c.stamp);
end
