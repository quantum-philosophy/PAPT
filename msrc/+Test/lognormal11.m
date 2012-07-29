init;

n_mu = 5;
n_sigma = 0.2;

fprintf('Y = exp(X)\n');
fprintf('X ~ N(%.2f, %.2f)\n\n', n_mu, n_sigma);

order = 3;
samples = 2e4;

fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);
fprintf('\n');

f = @(x) exp(n_mu + n_sigma * x);

%
% Exact solution
%

l_mu = exp(n_mu + n_sigma^2/2);
l_var = (exp(n_sigma^2) - 1) * exp(2 * n_mu + n_sigma^2);

fprintf('Exact stats:\n');
fprintf('mu = %.2f, var = %.2f\n\n', l_mu, l_var);

%
% Monte-Carlo
%

fprintf('Monte-Carlo simulation...');

t = tic;
[ mu, var, out_MC ] = MonteCarlo.sample1D(f, [ 1 1 ], samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PolynomialChaos([ 1 1 ], order);
fprintf(' %.2f seconds.\n', toc(t));

points = pc.points;
fprintf('Number of quadrature points: %d\n', points);

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out_PC ] = pc.sample(f, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Plots
%

Utils.compare(out_MC, out_PC, { 'MC', 'PC' }, ...
  @(x, i) lognpdf(x, n_mu, n_sigma));
