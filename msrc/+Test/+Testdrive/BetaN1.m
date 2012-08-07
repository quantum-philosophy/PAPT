init;

sdim = 4;

order = 5;
samples = 1e5;

alpha = 2;
beta = 2;
a = -1;
b = 1;

f = @(x) exp(prod(x, 1));

%% Monte-Carlo
%

out_MC = f(ibetarnd(alpha, beta, a, b, sdim, samples)).';

e = mean(out_MC);
v = var(out_MC);

fprintf('MC expectation: %.4f\n', e);
fprintf('MC variance: %.4f\n', v);

%% Polynomial Chaos
%

method = struct( ...
  'chaosOrder', order, ...
  'quadratureName', 'GaussJacobi', ...
  'quadratureType', 'Adaptive', ...
  'quadratureLevel', 5, ...
  'jacobiAlpha', alpha - 1, ...
  'jacobiBeta', beta - 1, ...
  'jacobiA', a, ...
  'jacobiB', b);

pc = PolynomialChaos.Jacobi([ sdim, 1 ], method);
[ e, v, out_PC ] = pc.sample(f, samples);

fprintf('PC expectation: %.4f\n', e);
fprintf('PC variance: %.4f\n', v);

%% Comparison
%

Utils.compareHistogram(out_MC, out_PC, ...
  { 'Monte Carlo', 'Polynomial Chaos' });
