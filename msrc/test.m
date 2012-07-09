init;

n_mu = 5;
n_sigma = 0.3;

f = @(x) exp(n_mu + n_sigma * x);

order = 3;
dimension = 1;
samples = 20000;

figure;

%
% Exact solution
%

l_mu = exp(n_mu + n_sigma^2/2);
l_var = (exp(n_sigma^2) - 1) * exp(2 * n_mu + n_sigma^2);

fprintf('Exact stats:\n');
fprintf('mu = %.2f, var = %.2f\n\n', l_mu, l_var);

x = linspace(0, 600, 100);
y = lognpdf(x, n_mu, n_sigma);
line(x, y, 'Color', 'r');

%
% Monte-Carlo
%

fprintf('Monte-Carlo simulation...');

t = tic;
[ mu, var, out ] = MC.perform(f, dimension, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

[ density, x ] = ksdensity(out);
line(x, density, 'Color', 'g');

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PC(dimension, order);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out ] = pc.perform(f, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

[ density, x ] = ksdensity(out);
line(x, density, 'Color', 'b');

legend('Exact', 'MC', 'PC');
