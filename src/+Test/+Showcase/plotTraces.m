function plotTraces
  setup;

  options = Test.configure('processorCount', 4);
  display(options);

  temperature = Temperature.Analytical.Transient(options);

  processors = 1:options.processorCount;

  Vnom = LeakagePower.Base.Vnom;
  scale = [ 1.1, 0.95, 0.90 ];

  Plot.figure(1000, 400);
  for i = 1:length(scale)
    T = Utils.toCelsius(temperature.compute(options.dynamicPower, ...
      'V', scale(i) * Vnom * ones(options.processorCount, 1)));
    for j = processors
      line(options.timeLine, T(j, :), 'Color', Color.pick(i), 'LineWidth', 1);
    end
  end
  Plot.label('Time, s', 'Temperature, C');
  Plot.limit(options.timeLine);
end
