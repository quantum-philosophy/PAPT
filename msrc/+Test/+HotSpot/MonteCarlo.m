function MonteCarlo
  clear all;
  setup;

  chaosSampleCount = 1e4;
  carloSampleCount = 1e3;

  options = configure;

  chaosOptions = Options('order', 6, ...
    'quadratureOptions', Options( ...
      'method', 'sparse', ...
      'ruleName', 'GaussHermiteHW', ...
      'order', 7));

  %
  % One polynomial chaos.
  %
  hotspot = HotSpot.Chaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, chaosOptions);

  display(hotspot);

  tic;
  [ Texp1, Tvar1, coefficients ] = hotspot.computeWithLeakage( ...
    options.powerProfile, options.leakage);
  fprintf('Polynomial chaos: %.2f s\n', toc);

  Tdata1 = hotspot.sample(coefficients, chaosSampleCount);

  %
  % Monte Carlo simulations.
  %
  hotspot = HotSpot.MonteCarlo(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, ...
    'sampleCount', carloSampleCount);

  display(hotspot);

  tic;
  [ Texp2, Tvar2, Tdata2 ] = hotspot.computeWithLeakageInParallel( ...
    options.powerProfile, options.leakage);
  fprintf('Monte Carlo: %.2f s\n', toc);

  Tdata2 = permute(Tdata2, [ 3 1 2 ]);

  time = 1e-3 * (0:(options.stepCount - 1));

  Utils.drawTemperature(time, ...
    { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
    { Tvar1, Tvar2 });

  k = floor(2 * options.stepCount / 3);

  Tdata1 = Utils.toCelsius(Tdata1(:, :, k));
  Tdata2 = Utils.toCelsius(Tdata2(:, :, k));

  Data.compare(Tdata2, Tdata1, ...
    'method', 'smooth', 'range', '4sigma', ...
    'layout', 'separate', 'draw', true, ...
    'labels', { 'Monte Carlo', 'Polynomial chaos' });
end
