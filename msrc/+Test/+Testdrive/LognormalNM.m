init;

sdim = 2;
ddim = 2;
order = 4;
samples = 2e4;

mu = [ 2, 2 ]';
sigma = [ 0.2, 0.8 ]';

A = [ 1, 0;
      0, 1 ];

fprintf('Y = exp((f1(X1, ..., X%d), ..., f%d(X1, ..., X%d))^T)\n', ...
  sdim, ddim, sdim);
for i = 1:sdim
  fprintf('X%d ~ N(%.2f, %.2f)\n', i, mu(i), sigma(i));
end
fprintf('\n');
fprintf('Number of spacial coordinates: %d\n', ddim);
fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);

%% Monte Carlo sampling
%
f = @(x) exp(A * (mu + sigma .* x));
[ ~, ~, mcRaw ] = MonteCarlo.sample1D(f, [ sdim ddim ], samples);

%% Polynomial Chaos expansion
%
pc = PolynomialChaos.Hermite([ sdim ddim ], struct('chaosOrder', order));

points = pc.points;

mu = irep(mu, 1, points);
sigma = irep(sigma, 1, points);

f = @(x) exp(A * (mu + sigma .* x));

[ ~, ~, pcRaw ] = pc.sample(f, samples);

%% Comparison
%
Utils.compare(mcRaw, pcRaw, 'draw', true);
