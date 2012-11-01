setup;

sampleCount = 1e5;
processorCount = [ 2 4 8 16 32 ];

for i = 1:length(processorCount)
  floorplan = File.join( ...
    File.trace, '..', 'Assets', sprintf('%02d.flp', processorCount(i)));

  PCA = ProcessVariation.analyze(floorplan, LeakagePower.Lnom);

  [ ~, dimension ] = size(PCA);

  fprintf('Processors: %d, Random variables: %d\n', ...
    processorCount(i), dimension);

  fprintf('%10s%20s%15s\n', 'Processor', 'Expectation, nm', 'Deviation, nm');

  L = PCA * normrnd(0, 1, dimension, sampleCount);
  expectation = mean(L, 2) * 1e9;
  deviation = sqrt(var(L, 0, 2)) * 1e9;

  for j = 1:processorCount(i)
    fprintf('%10d%20.2f%15.2f\n', j, expectation(j), deviation(j));
  end

  fprintf('\n');
end