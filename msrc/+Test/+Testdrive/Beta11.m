init;

order = 1;
samples = 1e5;

alpha = 1;
beta = 1;
a = 0;
b = 1;

e = (alpha * b + beta * a) / (alpha + beta);
v = (alpha * beta * (b - a)^2) / (alpha + beta)^2 / (alpha + beta + 1);

fprintf('Exact expectation: %.4f\n', e);
fprintf('Exact variance: %.4f\n', v);

%% Monte-Carlo
%

% out_MC = ibetarnd(alpha, beta, a, b, 1, samples).';
out_MC = betarnd(alpha, beta, 1, samples).';

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

pc = PolynomialChaos.Jacobi([ 1 1 ], order, method);
[ e, v, out_PC ] = pc.sample(@(x) x, samples);

fprintf('PC expectation: %.4f\n', e);
fprintf('PC variance: %.4f\n', v);

%% Comparison
%

Utils.compareSmooth(out_MC, out_PC, struct(...
  'labels', {{ 'MC', 'PC', 'Exact' }}, 'x', (a - 1):0.01:(b + 1), ...
  'exact', @(x, i) ibetapdf(x, alpha, beta, a, b)));
