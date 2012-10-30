function Chaos
  clear all;
  setup;

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
  [ Texp, Tvar ] = hotspot.computeWithLeakage(Pdyn, leakage);
  fprintf('Polynomial chaos expansion: %.2f s\n', toc);

  time = 1e-3 * (1:size(Pdyn, 2));

  draw(time, Utils.toCelsius(Texp), Tvar);

  Plot.title('Temperature');
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(time);
end

function draw(time, expectation, variance)
  figure;

  deviation = sqrt(variance);
  processorCount = size(expectation, 1);

  for i = 1:processorCount
    color = Color.pick(i);
    line(time, expectation(i, :), ...
      'Color', color, 'LineWidth', 1.5);
    line(time, expectation(i, :) + 3 * deviation(i, :), ...
      'Color', color, 'LineStyle', '--');
  end
end
