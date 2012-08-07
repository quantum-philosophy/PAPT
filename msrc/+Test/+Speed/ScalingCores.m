init;

steps = 1000;

%% Define the needed measurements.
%
X = [ 2 4 8 16 32 ];

fprintf('%15s%15s%15s%15s\n', 'Cores', 'Chaos, s', 'Kutta, h', 'Speedup, x');

Y = zeros(length(X), 2);

for i = 1:length(X)
  c = Config('cores', X(i), 'constructionMethod.chaosOrder', 4, 'steps', steps);

  fprintf('%15d', c.cores);

  %% Initialize the solver.
  %
  ch = HotSpot.Chaos(c.hotspotArguments{:}, c.chaosArguments{:});
  kt = HotSpot.Kutta(c.hotspotArguments{:});

  ch.configureLeakage(c.dynamicPower);
  kt.configureLeakage(c.dynamicPower);

  %% Perform the analysis.
  %
  t = tic;
  ch.solve(c.dynamicPower);
  Y(i, 1) = toc(t);

  fprintf('%15.2f', Y(i, 1));

  t = tic;
  kt.solve(c.dynamicPower, zeros(kt.sdim, 1));
  Y(i, 2) = toc(t) * c.monteCarloSamples;

  fprintf('%15.2f', Y(i, 2) / 60 / 60);

  fprintf('%15.2e', Y(i, 2) / Y(i, 1));
  fprintf('\n');
end

figure;
line(X, Y(:, 1), 'Color', Utils.pickColor(1), 'Marker', 'o')
line(X, Y(:, 2), 'Color', Utils.pickColor(2), 'Marker', 'x')
title('Scaling with Number of Cores');
xlabel('Steps');
ylabel('Time, s');

legend('Polynomial Chaos', 'Monte Carlo');
