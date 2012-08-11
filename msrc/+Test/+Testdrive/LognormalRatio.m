init;

sdim = 2;
order = 4;
samples = 1e5;

mu = [ 1, 1 ]';
sigma = [ 0.5, 0.5 ]';

fprintf('Z = X / Y\n');
fprintf('X ~ LogN(%.2f, %.2f)\n', mu(1), sigma(1));
fprintf('Y ~ LogN(%.2f, %.2f)\n', mu(2), sigma(2));
fprintf('\n');
fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);

%% Monte Carlo sampling
%
f = @(x) exp(mu(1) + sigma(1) * x(1))./exp(mu(2) + sigma(2) * x(2));

[ ~, ~, mcRaw ] = MonteCarlo.sample1D(f, [ sdim 1 ], samples);

%% Polynomial Chaos expansion
%
pc = PolynomialChaos.Hermite([ sdim 1 ], struct('chaosOrder', order));

mu = irep(mu, 1, pc.points);
sigma = irep(sigma, 1, pc.points);

f = @(x) exp(mu(1, :) + sigma(1, :) .* x(1, :)) ./ ...
  exp(mu(2, :) + sigma(2, :) .* x(2, :));

[ ~, ~, pcRaw ] = pc.sample(f, samples);

%% Comparison
%
Utils.compare(mcRaw, pcRaw, ...
  'method', 'smooth', 'range', '3sigma', ...
  'function', 'pdf', 'draw', true);
