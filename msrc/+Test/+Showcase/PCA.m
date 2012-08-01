init;

samples = 10^5;
coreSet = [ 2, 4, 8, 16, 32 ];
coreCount = length(coreSet);

for i = 1:coreCount
  floorplan = Utils.resolveTest(coreSet(i));

  P = PrincipalComponent.perform(floorplan);

  [ cores, sdim ] = size(P);

  assert(cores == coreSet(i), 'The constructed matrix is invalid.');

  fprintf('Cores: %d, Random variables: %d\n', cores, sdim);
  fprintf('%10s%20s%15s\n', 'Core', 'Expectation, nm', 'Deviation, nm');

  L = P * normrnd(0, 1, sdim, samples);
  exp = mean(L, 2) * 1e9;
  std = sqrt(var(L, 0, 2)) * 1e9;

  for j = 1:cores
    fprintf('%10d%20.2f%15.2f\n', j, exp(j), std(j));
  end

  fprintf('\n');
end
