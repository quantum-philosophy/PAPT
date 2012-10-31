function MonteCarlo
  clear all;
  setup;

  options = Test.HotSpot.configure;

  chaosOptions = Options('order', 10, ...
    'quadratureOptions', Options('order', 11));

  %
  % One polynomial chaos.
  %
  hotspot = HotSpot.StepwiseChaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, chaosOptions);

  display(hotspot);

  tic;
  [ Texp1, Tvar1 ] = hotspot.computeWithLeakage( ...
    options.powerProfile, options.leakage);
  fprintf('Polynomial chaos: %.2f s\n', toc);

  %
  % Monte Carlo simulations.
  %
  hotspot = HotSpot.MonteCarlo(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine);

  display(hotspot);

  tic;
  [ Texp2, Tvar2 ] = hotspot.computeWithLeakage( ...
    options.powerProfile, options.leakage);
  fprintf('Monte Carlo: %.2f s\n', toc);

  time = 1e-3 * (1:options.stepCount);

  Test.HotSpot.draw(time, ...
    { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
    { Tvar1, Tvar2 });
end
