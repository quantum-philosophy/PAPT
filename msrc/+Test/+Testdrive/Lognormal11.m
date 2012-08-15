init;

order = 3;
samples = 2e4;

mu = 5;
sigma = 0.2;

fprintf('Y = exp(X)\n');
fprintf('X ~ N(%.2f, %.2f)\n', mu, sigma);
fprintf('\n');
fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);

f = @(x) exp(mu + sigma * x);

%% Monte Carlo sampling
%
[ ~, ~, mcRaw ] = MonteCarlo.sample1D(f, [ 1 1 ], samples);

%% Polynomial Chaos expansion
%
pc = PolynomialChaos.Hermite([ 1 1 ], struct('chaosOrder', order));
[ ~, ~, pcRaw ] = pc.sample(f, samples);

%% Comparison
%
Stats.compare(mcRaw, pcRaw, 'draw', true);
