init;

X = [ 2, 4, 8, 16, 32 ];

fprintf('%15s%15s\n', 'Cores', 'PC + 1');

for i = 1:length(X)
  cores = X(i);

  floorplan = Utils.resolveTest(cores);

  C = PrincipalComponent.perform(floorplan);

  assert(size(C, 1) == cores, 'The constructed matrix is invalid.');

  fprintf('%15d%15d\n', cores, size(C, 2));
end
