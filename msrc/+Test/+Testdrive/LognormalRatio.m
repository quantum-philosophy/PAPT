init;

sdim = 2;

n_mu = [ 2, 2 ]';
n_sigma = [ 0.7, 0.7 ]';

fprintf('Z = X / Y\n');
fprintf('X ~ LogN(%.2f, %.2f)\n', n_mu(1), n_sigma(1));
fprintf('Y ~ LogN(%.2f, %.2f)\n', n_mu(2), n_sigma(2));
fprintf('\n');

order = 10;
samples = 1e5;

fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);
fprintf('\n');

%
% Monte-Carlo
%

fprintf('Monte-Carlo simulation...');

f = @(x) exp(n_mu(1) + n_sigma(1) * x(1))./exp(n_mu(2) + n_sigma(2) * x(2));

t = tic;
[ mu, var, out_MC ] = MonteCarlo.sample1D(f, [ sdim 1 ], samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PolynomialChaos.Hermite([ sdim 1 ], struct('chaosOrder', order));
fprintf(' %.2f seconds.\n', toc(t));

points = pc.points;
fprintf('Number of quadrature points: %d\n', points);

n_mu = irep(n_mu, 1, points);
n_sigma = irep(n_sigma, 1, points);

f = @(x) exp(n_mu(1, :) + n_sigma(1, :) .* x(1, :))./exp(n_mu(2, :) + n_sigma(2, :) .* x(2, :));

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out_PC ] = pc.sample(f, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Plots
%

Utils.compareSmooth(...
  out_MC, out_PC, struct('labels', {{ 'MC', 'PC' }}));
