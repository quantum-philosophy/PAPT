function detectThermalRunaway
  close all;
  setup;

  options = Test.configure('processorCount', 2);

  chaos = Temperature.Chaos.DynamicSteadyState( ...
    options.temperatureOptions, options.chaosOptions);
  display(chaos);

  dimensionCount = chaos.process.dimensionCount;

  pairs = combnk(1:dimensionCount, 2);

  x = linspace(-5, 5, 20);
  y = x;

  [ X, Y ] = meshgrid(x, y);

  for k = 1:size(pairs, 1)
    figure;
    Plot.title('Variables %d and %d', pairs(k, 1), pairs(k, 2));

    rvs = zeros(dimensionCount, length(x)^2);
    rvs(pairs(k, :), :) = [ X(:)'; Y(:)' ];

    [ T, output ] = chaos.solve(options.dynamicPower, rvs, ...
      'iterationLimit', 20, 'temperatureLimit', Utils.toKelvin(200), ...
      'tolerance', 0.5);

    R = reshape(output.iterationCount, length(x), []);

    surfc(X, Y, R);
    view(0, 90);
  end
end
