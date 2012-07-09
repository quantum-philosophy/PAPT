init;

n_mu = 5;
n_sigma = 0.3;

f = @(x) exp(n_mu + n_sigma * x);

l_mu = exp(n_mu + n_sigma^2/2);
l_sigma2 = (exp(n_sigma^2) - 1) * exp(2 * n_mu + n_sigma^2);

dimension = 1;

fprintf('Exact stats:\n');
fprintf('mu = %.2f, sigma2 = %.2f\n\n', l_mu, l_sigma2);

fprintf('Monte-Carlo simulation...');

tic;
[ mu, sigma2, out ] = MC.perform(f, dimension);
fprintf(' %.2f seconds.\n', toc);

fprintf('mu = %.2f, sigma2 = %.2f\n\n', mu, sigma2);

% hist(out, 50);

fprintf('Polynomial Chaos preparation...');
tic
pc = PC(dimension);
fprintf(' %.2f seconds.\n', toc);

fprintf('Polynomial Chaos simulation...');
tic
[ mu, sigma2 ] = pc.perform(f);
fprintf(' %.2f seconds.\n', toc);

fprintf('mu = %.2f, sigma2 = %.2f\n\n', mu, sigma2);
