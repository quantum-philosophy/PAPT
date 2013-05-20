setup;
close all;

sampleCount = 1e5;
processorCount = [ 2 4 8 16 32 ];

for i = 1:length(processorCount)
  options = Test.configure('processorCount', processorCount(i));

  process = ProcessVariation.Normal(options.processOptions);

  title = sprintf('Processors: %d, Variables: %d\n', ...
    processorCount(i), process.dimensionCount);

  plot(options.die);
  Plot.title(title);

  fprintf(title);
  fprintf('%10s%20s%15s\n', 'Processor', 'Expectation, nm', 'Deviation, %');

  L = process.sample(sampleCount);
  expectation = mean(L, 1) * 1e9;
  deviation = sqrt(var(L, 0, 1)) * 1e9;

  for j = 1:processorCount(i)
    fprintf('%10d%20.2f%15.2f\n', j, expectation(j), ...
      deviation(j) / expectation(j) * 100);
  end

  fprintf('\n');
end
