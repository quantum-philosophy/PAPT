function Chaos
  clear all;
  setup;

  use('Vendor', 'heatmaps');

  options = configure;

  %
  % One polynomial chaos.
  %
  hotspot = HotSpot.Chaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, options.chaosOptions);

  display(hotspot);

  tic;
  [ Texp1, Tvar1, coefficients1 ] = ...
    hotspot.computeWithLeakage(options.powerProfile, options.leakage);
  fprintf('Polynomial chaos: %.2f s\n', toc);

  %
  % Another polynomial chaos.
  %
  hotspot = HotSpot.StepwiseChaos(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine, options.chaosOptions);

  tic;
  [ Texp2, Tvar2, coefficients2 ] = ...
    hotspot.computeWithLeakage(options.powerProfile, options.leakage);
  fprintf('Stepwise polynomial chaos: %.2f s\n', toc);

  time = 1e-3 * (1:options.stepCount);

  Utils.drawTemperature(time, ...
    { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
    { Tvar1, Tvar2 }, 'labels', { 'PC 1', 'PC stepwise' });

  showCoefficients(time, { coefficients1, coefficients2 });
end

function showCoefficients(~, coefficientSet)
  setCount = length(coefficientSet);
  [ ~, processorCount, ~ ] = size(coefficientSet{1});

  for i = 1:processorCount
    figure;
    for j = 1:setCount
      subplot(1, setCount, j);
      heatmap(flipud(abs(squeeze(coefficientSet{j}(2:end, i, :)))));
      Plot.title('Magnitude (PE%d)', i);
      Plot.label('Time', 'Coefficient');
    end
  end
end
