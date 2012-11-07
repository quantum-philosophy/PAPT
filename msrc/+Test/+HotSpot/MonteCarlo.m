function MonteCarlo
  clear all;
  setup;

  chaosSampleCount = 1e5;
  carloSampleCount = 1e5;

  options = configure;

  chaosOptions = Options('order', 6, ...
    'quadratureOptions', Options( ...
      'method', 'sparse', ...
      'ruleName', 'GaussHermiteHW', ...
      'order', 6 + 1));

  time = options.samplingInterval * (0:(options.stepCount - 1));

  timeSlice = 0.04;
  k = floor(timeSlice / options.samplingInterval);

  %
  % One polynomial chaos.
  %
  chaos = HotSpot.Chaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, chaosOptions);

  display(chaos);

  %
  % Monte Carlo simulations.
  %
  mc = HotSpot.MonteCarlo(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, ...
    'sampleCount', carloSampleCount, 'verbose', 'true');

  display(mc);

  tic;
  [ Texp1, Tvar1, coefficients ] = chaos.computeWithLeakage( ...
    options.powerProfile, options.leakage);
  fprintf('Polynomial chaos: construction time %.2f s.\n', toc);

  %
  % Comparison of expectations, variances, and PDFs.
  %
  if Terminal.question('Compare expectations, variances, and PDFs? ')
    tic;
    Tdata1 = chaos.sample(coefficients, chaosSampleCount);
    fprintf('Polynomial chaos: sampling time %.2f s (%d samples).\n', ...
      toc, chaosSampleCount);

    [ Texp2, Tvar2, Tdata2 ] = mc.computeWithLeakageInParallel( ...
      options.powerProfile, options.leakage);

    labels = { 'PC', 'MC' };

    Utils.drawTemperature(time, ...
      { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
      { Tvar1, Tvar2 }, 'labels', labels);

    Tdata1 = Utils.toCelsius(Tdata1(:, :, k));
    Tdata2 = Utils.toCelsius(Tdata2(:, :, k));

    Data.compare(Tdata1, Tdata2, ...
      'method', 'histogram', 'range', 'unbounded', ...
      'layout', 'separate', 'draw', true, ...
      'labels', labels);
  end

  %
  % Sweeping of the random parameters.
  %
  rvs = -7:0.2:7;
  index = uint8(1);
  while Terminal.question('Sweep random variables? ')
    index = Terminal.request( ...
     'prompt', sprintf('Which random variables? [%s] ', Utils.toString(index)), ...
     'type', 'uint8', 'default', index);

     if any(index > chaos.dimension), continue; end

    RVs = zeros(length(rvs), options.processorCount + 1);
    for i = index
      RVs(:, i) = rvs;
    end

    Tdata1 = chaos.evaluate(coefficients, RVs);

    Tdata2 = mc.evaluateWithLeakageInParallel( ...
      options.powerProfile, options.leakage, RVs);

    Tdata1 = Utils.toCelsius(Tdata1(:, :, k));
    Tdata2 = Utils.toCelsius(Tdata2(:, :, k));

    figure;
    for i = 1:options.processorCount
      color = Color.pick(i);
      line(rvs, Tdata1(:, i), 'Color', color, 'Marker', 'o');
      line(rvs, Tdata2(:, i), 'Color', color, 'Marker', 'x');
    end
    Plot.title('Sweep at %.3f s', timeSlice);
    Plot.label('Random parameter', 'Temperature, C');
  end
end
