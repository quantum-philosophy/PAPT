init;

sdim = 2;
order = 4;
samples = 2e4;

mu = [ 2, 2 ]';
sigma = [ 0.2, 0.8 ]';

fprintf('Y = exp(X1 + ... + X%d)\n', sdim);
for i = 1:sdim
  fprintf('X%d ~ N(%.2f, %.2f)\n', i, mu(i), sigma(i));
end
fprintf('\n');
fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);

%% Monte Carlo sampling
%
f = @(x) exp(sum(mu + sigma .* x));

[ ~, ~, mcRaw ] = MonteCarlo.sample1D(f, [ sdim 1 ], samples);

%% Polynomial Chaos expansion
%
pc = PolynomialChaos.Hermite([ sdim 1 ], struct('chaosOrder', order));

points = pc.points;

mu = irep(mu, 1, points);
sigma = irep(sigma, 1, points);

f = @(x) exp(sum(mu + sigma .* x));

[ ~, ~, pcRaw ] = pc.sample(f, samples);

%% Comparison
%
Utils.compare(mcRaw, pcRaw, 'draw', true);
