function compareSpeed(varargin)
  setup;

  sampleCount = 1e4;

  experiment = Options(varargin{:});
  experimentCount = length(experiment.range);
  repeat = experiment.get('repeat', ones(1, experimentCount));

  fprintf('%15s%15s%15s%15s%15s%15s\n', ...
    experiment.shortName, 'Chaos, s', ...
    'Analytic, m', 'Speedup, x', ...
    'Numeric, h', 'Speedup, x');

  measurements = zeros(experimentCount, 3);

  for i = 1:experimentCount
    parameter = experiment.range(i);

    fprintf('%15d', parameter);

    options = experiment.configure(parameter);

    chaos = Temperature.Chaos.Transient(options);
    analytic = Temperature.Analytical.Transient(options.temperatureOptions);
    numeric = Temperature.Numerical.Transient(options.temperatureOptions);

    for j = 1:repeat(i)
      tic;
      analytic.compute(options.dynamicPower);
      measurements(i, 2) = measurements(i, 2) + toc * sampleCount;

      tic;
      numeric.compute(options.dynamicPower);
      measurements(i, 3) = measurements(i, 3) + toc * sampleCount;

      tic;
      chaos.compute(options.dynamicPower);
      measurements(i, 1) = measurements(i, 1) + toc;
    end

    measurements(i, :) = measurements(i, :) / repeat(i);

    fprintf('%15.2f', measurements(i, 1));
    fprintf('%15.2f', measurements(i, 2) / 60);
    fprintf('%15.2e', measurements(i, 2) / measurements(i, 1));
    fprintf('%15.2f', measurements(i, 3) / 60 / 60);
    fprintf('%15.2e', measurements(i, 3) / measurements(i, 1));
    fprintf('\n');
  end

  figure;

  line(experiment.range, measurements(:, 1), ...
    'Color', Color.pick(1), 'Marker', 'o')
  line(experiment.range, measurements(:, 2), ...
    'Color', Color.pick(2), 'Marker', 'x')
  line(experiment.range, measurements(:, 3), ...
    'Color', Color.pick(3), 'Marker', 's')

  set(gca, 'YScale', 'log');

  Plot.title('Comparison of computational speed: %s', experiment.name);
  Plot.label(experiment.name, 'log(Time, s)');
  Plot.legend('Polynomial Chaos', ...
    'Monte Carlo (Analytic)', 'Monte Carlo (Numeric)');
end
