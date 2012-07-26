init;

sdim = 2;
ddim = 2;

n_mu = [ 2, 2 ]';
n_sigma = [ 0.2, 0.8 ]';

A = [ 1, 0;
      0, 1 ];

fprintf('Y = exp((f1(X1, ..., X%d), ..., f%d(X1, ..., X%d))^T)\n', ...
  sdim, ddim, sdim);
for i = 1:sdim
  fprintf('X%d ~ N(%.2f, %.2f)\n', i, n_mu(i), n_sigma(i));
end
fprintf('\n');

order = 4;
samples = 2e4;

fprintf('Number of spacial coordinates: %d\n', ddim);
fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);
fprintf('\n');

%
% Monte-Carlo
%

fprintf('Monte-Carlo simulation...');

f = @(x) exp(A * (n_mu + n_sigma .* x));

t = tic;
[ mu, var, out_MC ] = MonteCarlo.perform(f, [ sdim ddim ], samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('%10s%10s\n', 'mu', 'var');
for i = 1:ddim
  fprintf('%10.2f%10.2f\n', mu(i), var(i, i));
end
fprintf('\n');

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PolynomialChaos([ sdim ddim ], order);
fprintf(' %.2f seconds.\n', toc(t));

points = pc.gq.points;
fprintf('Number of quadrature points: %d\n', points);

n_mu = irep(n_mu, 1, points);
n_sigma = irep(n_sigma, 1, points);

f = @(x) exp(A * (n_mu + n_sigma .* x));

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out_PC ] = pc.perform(f, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('%10s%10s\n', 'mu', 'var');
for i = 1:ddim
  fprintf('%10.2f%10.2f\n', mu(i), var(i, i));
end
fprintf('\n');

%
% Plots
%

Utils.compare(out_MC, out_PC, { 'MC', 'PC' });
