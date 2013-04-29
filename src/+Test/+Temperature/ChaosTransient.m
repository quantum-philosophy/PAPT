function ChaosTransient
  close all;
  setup;

  options = Test.configure;

  plot(options.schedule);
  Utils.plotFloorplan(options.temperatureOptions.floorplan);

  chaos = Temperature.Chaos.Transient( ...
    options.temperatureOptions, options.chaosOptions);

  display(chaos);

  tic;
  [ Texp, output ] = chaos.compute(options.dynamicProfile);
  fprintf('Polynomial chaos: %.2f s\n', toc);

  time = options.samplingInterval * (1:options.stepCount);

  Utils.drawTemperature(time, { Utils.toCelsius(Texp) }, { output.Tvar });
  showCoefficients(time, { output.coefficients });
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
