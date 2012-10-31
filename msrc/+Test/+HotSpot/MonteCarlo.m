function MonteCarlo
  clear all;
  setup;

  [ hotspotArguments, Pdyn, leakage ] = Test.HotSpot.configure;

  chaosOptions = Options('order', 10, ...
    'quadratureOptions', Options('order', 11));

  %
  % One polynomial chaos.
  %
  hotspot = HotSpot.Chaos(hotspotArguments{:}, chaosOptions);

  display(hotspot);

  tic;
  [ Texp1, Tvar1 ] = hotspot.computeWithLeakage(Pdyn, leakage);
  fprintf('Polynomial chaos: %.2f s\n', toc);

  %
  % Monte Carlo simulations.
  %
  hotspot = HotSpot.MonteCarlo(hotspotArguments{:});

  display(hotspot);

  tic;
  [ Texp2, Tvar2 ] = hotspot.computeWithLeakage(Pdyn, leakage);
  fprintf('Monte Carlo: %.2f s\n', toc);

  time = 1e-3 * (1:size(Pdyn, 2));

  Test.HotSpot.draw(time, ...
    { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
    { Tvar1, Tvar2 });
end
