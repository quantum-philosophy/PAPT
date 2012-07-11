init;

sdim = 2;

n_mu = [ 2, 2 ]';
n_sigma = [ 0.2, 0.8 ]';

fprintf('Y = exp(X1 + ... + X%d)\n', sdim);
for i = 1:sdim
  fprintf('X%d ~ N(%.2f, %.2f)\n', i, n_mu(i), n_sigma(i));
end
fprintf('\n');

order = 4;
samples = 2e4;

fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);
fprintf('\n');

%
% Monte-Carlo
%

fprintf('Monte-Carlo simulation...');

f = @(x) exp(sum(n_mu + n_sigma .* x));

t = tic;
[ mu, var, out_MC ] = MonteCarlo.perform(f, [ sdim 1 ], samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PolynomialChaos(sdim, order);
fprintf(' %.2f seconds.\n', toc(t));

count = pc.cq.count;

n_mu = irep(n_mu, 1, count);
n_sigma = irep(n_sigma, 1, count);

f = @(x) exp(sum(n_mu + n_sigma .* x));

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out_PC ] = pc.perform(f, [ sdim 1 ], samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Plots
%

Utils.compare(out_MC, out_PC, { 'MC', 'PC' });
