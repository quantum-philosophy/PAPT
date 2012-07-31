init;

c = Config('polynomialOrder', 4, 'cores', 4);

%% Initialize the solver.
%
ch = HotSpot.Chaos(c.hotspotArguments{:}, c.polynomialOrder, c.integrationMethod);
kt = HotSpot.Kutta(c.hotspotArguments{:});

%% Define the needed measurements.
%
X = [ 10, 100, 1000, 10000, 100000 ];

fprintf('%15s%15s%15s%15s\n', 'Steps', 'Chaos, s', 'Kutta, h', 'Speedup, x');

Y = zeros(length(X), 2);

for i = 1:length(X)
  c.tune('steps', X(i));

  fprintf('%15d', c.steps);

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
title('Scaling with Number of Steps');
xlabel('Steps');
ylabel('Time, s');

legend('Polynomial Chaos', 'Monte Carlo');
