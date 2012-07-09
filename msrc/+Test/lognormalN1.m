init;

dimension = 2;

n_mu = [ 2, 2 ];
n_sigma = [ 0.2, 0.8 ];

fprintf('Y = exp(X1 + ... + X%d)\n', dimension);
for i = 1:dimension
  fprintf('X%d ~ N(%.2f, %.2f)\n', i, n_mu(i), n_sigma(i));
end
fprintf('\n');

order = 4;
samples = 2e4;

fprintf('Order of PC: %d\n', order);
fprintf('Number of samples: %d\n', samples);
fprintf('\n');

f = @(x) exp(sum(n_mu + n_sigma .* x));

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

density_MC = ksdensity(out_MC, x);
line(x, density_MC, 'Color', 'b');

density_PC = ksdensity(out_PC, x);
line(x, density_PC, 'Color', 'g');

title('PDF');
legend('MC', 'PC');
xlim([ x(1), x(end) ]);

subplot(2, 1, 2);

line(x, density_MC - density_PC, 'Color', 'r');

title('Error');
legend('MC - PC');
xlim([ x(1), x(end) ]);

%
% RMSE
%

n = length(x);

rmse = sqrt(sum((density_MC - density_PC).^2) / n);

fprintf('RMSE: %.4e\n', rmse);
