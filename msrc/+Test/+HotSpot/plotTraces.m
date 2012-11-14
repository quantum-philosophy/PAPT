function motivation2
  options = configure;

  mc = HotSpot.MonteCarlo(options.floorplan, ...
    options.hotspotConfig, options.hotspotLine);

  processors = [ 1, 2 ];

  Lnom = LeakagePower.Lnom;
  scale = [ 1, 0.90, 0.85 ];

  t = (0:(options.stepCount - 1)) * options.samplingInterval;

  figure;
  for i = 1:length(scale)
    T = Utils.toCelsius(mc.computeWithLeakage( ...
      options.powerProfile, options.leakage, scale(i) * Lnom));
    for j = processors
      line(t, T(j, :), 'Color', Color.pick(i), 'LineWidth', 1);
    end
  end
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(t);
end
