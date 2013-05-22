function plotTraces
  setup;

  options = Test.configure('processorCount', 4, 'powerScale', 1);

  temperature = Temperature.Analytical.Transient(options);

  processors = [ 1 ];

  Lnom = LeakagePower.Base.Lnom;
  scale = [ 1, 0.95, 0.90 ];

  t = (0:(options.stepCount - 1)) * options.samplingInterval;

  figure;
  for i = 1:length(scale)
    T = Utils.toCelsius(temperature.compute( ...
      options.dynamicPower, 'L', scale(i) * Lnom));
    for j = processors
      line(t, T(j, :), 'Color', Color.pick(i), 'LineWidth', 1);
    end
  end
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(t);
end
