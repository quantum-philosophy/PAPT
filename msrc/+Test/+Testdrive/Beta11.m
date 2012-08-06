init;

order = 1;
samples = 1e5;

alpha = 2;
beta = 2;
a = -1;
b = 1;

e = (alpha * b + beta * a) / (alpha + beta);
v = (alpha * beta * (b - a)^2) / (alpha + beta)^2 / (alpha + beta + 1);

fprintf('Exact expectation: %.4f\n', e);
fprintf('Exact variance: %.4f\n', v);

f = @(x) x;

%% Monte-Carlo
%

out_MC = f(ibetarnd(alpha, beta, a, b, 1, samples)).';

e = mean(out_MC);
v = var(out_MC);

fprintf('MC expectation: %.4f\n', e);
fprintf('MC variance: %.4f\n', v);

%% Polynomial Chaos
%

method = struct( ...
  'quadratureName', 'GaussJacobi', ...
  'quadratureType', 'Tensor', ...
  'jacobiAlpha', alpha - 1, ...
  'jacobiBeta', beta - 1, ...
  'jacobiA', a, ...
  'jacobiB', b);

pc = PolynomialChaos.Jacobi([ 1, 1 ], order, method);
[ e, v, out_PC ] = pc.sample(f, samples);

fprintf('PC expectation: %.4f\n', e);
fprintf('PC variance: %.4f\n', v);

%% Comparison
%

Utils.compareHistogram(out_MC, out_PC)

x = xlim(gca);
x = linspace(x(1), x(2), 200);
h = line(x, ibetapdf(x, alpha, beta, a, b), ...
  'Color', Utils.pickColor(5), 'LineStyle', '--');
set(h,'erasemode','xor');
set(h,'erasemode','background');

legend('Monte Carlo', 'Polynomial Chaos', 'Exact');
