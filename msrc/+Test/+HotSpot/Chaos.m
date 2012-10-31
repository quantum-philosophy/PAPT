function Chaos
  clear all;
  setup;

  use('Vendor', 'heatmaps');

  path = File.join(File.trace, '..', 'Assets');

  floorplan = File.join(path, '04.flp');
  powerProfile = File.join(path, '04.ptrace');

  hotspotConfig = File.join(path, 'hotspot.config');
  hotspotLine = 'sampling_intvl 1e-3';

  leakageFilename = File.join(path, 'inverter_45nm.leak');

  Pdyn = 2 * dlmread(powerProfile, '', 1, 0).';

  leakage = LeakagePower(Pdyn, ...
    'filename', leakageFilename, ...
    'order', [ 1, 2 ], ...
    'scale', [ 1, 0.7, 0; 1, 1, 1 ]);

  hotspot = HotSpot.Chaos(floorplan, hotspotConfig, hotspotLine);

  display(hotspot);

  tic;
  [ Texp1, Tvar1, coefficients1 ] = ...
    hotspot.computeWithLeakage(Pdyn, leakage);
  fprintf('Polynomial chaos expansion: %.2f s\n', toc);

  tic;
  [ Texp2, Tvar2, coefficients2 ] = ...
    hotspot.computeWithLeakageStepwise(Pdyn, leakage);
  fprintf('Polynomial chaos expansion stepwise: %.2f s\n', toc);

  time = 1e-3 * (1:size(Pdyn, 2));

  draw(time, ...
    { Utils.toCelsius(Texp1), Utils.toCelsius(Texp2) }, ...
    { Tvar1, Tvar2 });

  showCoefficients(time, { coefficients1, coefficients2 });
end

function draw(time, expectationSet, varianceSet)
  setCount = length(expectationSet);
  processorCount = size(expectationSet{1}, 1);

  for i = 1:processorCount
    figure;
    for j = 1:setCount
      color = Color.pick(j);
      line(time, expectationSet{j}(i, :), ...
        'Color', color, 'LineWidth', 1);
      line(time, expectationSet{j}(i, :) + sqrt(varianceSet{j}(i, :)), ...
        'Color', color, 'LineStyle', '--');
    end
    Plot.title('Temperature (PE%d)', i);
    Plot.label('Time, s', 'Temperature, C');
    Plot.limit(time);
  end
end

function showCoefficients(time, coefficientSet);
  setCount = length(coefficientSet);
  [ coefficientCount, processorCount, stepCount ] = size(coefficientSet{1});

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
