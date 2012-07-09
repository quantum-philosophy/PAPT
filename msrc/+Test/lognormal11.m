init;

dimension = 1;

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
[ mu, var, out_MC ] = MC.perform(f, dimension, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Polynomial Chaos expansion
%

fprintf('Polynomial Chaos preparation...');
t = tic;
pc = PC(dimension, order);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('Polynomial Chaos simulation...');
t = tic;
[ mu, var, out_PC ] = pc.perform(f, samples);
fprintf(' %.2f seconds.\n', toc(t));

fprintf('mu = %.2f, var = %.2f\n\n', mu, var);

%
% Plots
%

min_MC = min(min(out_MC));
max_MC = max(max(out_MC));

min_PC = min(min(out_PC));
max_PC = max(max(out_PC));

x = linspace(max(min_MC, min_PC), min(max_MC, max_PC), 200);

figure;

subplot(2, 1, 1);

density = lognpdf(x, n_mu, n_sigma);
line(x, density, 'Color', 'r');

density_MC = ksdensity(out_MC, x);
line(x, density_MC, 'Color', 'b');

density_PC = ksdensity(out_PC, x);
line(x, density_PC, 'Color', 'g');

title('PDF');
legend('Exact', 'MC', 'PC');
xlim([ x(1), x(end) ]);

subplot(2, 1, 2);

line(x, density - density_MC, 'Color', 'b');
line(x, density - density_PC, 'Color', 'g');

title('Error');
legend('MC', 'PC');
xlim([ x(1), x(end) ]);

%
% RMSE
%

n = length(x);

rmse_MC = sqrt(sum((density_MC - density).^2) / n);
rmse_PC = sqrt(sum((density_PC - density).^2) / n);

fprintf('RMSE for MC: %.4e\n', rmse_MC);
fprintf('RMSE for PC: %.4e\n', rmse_PC);
