init;

monteCarloSamples = 1e4;
steps = 1000;

%% Define the needed measurements.
%
X = [ 2, 4, 8, 16, 32 ];

fprintf('%15s%15s%15s%15s\n', 'Cores', 'Chaos, s', 'Kutta, h', 'Speedup, x');

Y = zeros(length(X), 2);

for i = 1:length(X)
  cores = X(i);

  fprintf('%15d', cores);

  [ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(cores);

  %% Initialize the solver.
  %
  ch = HotSpot.Chaos(floorplan, config, configLine);
  kt = HotSpot.Kutta(floorplan, config, configLine);

  %% Read the specified template and fit it to the desired length.
  %
  P = dlmread(powerTrace, '', 1, 0)';
  PP = Utils.replicate(P, steps);

  %% Perform the analysis.
  %
  time = tic;
  ch.solve(PP);
  Y(i, 1) = toc(time);

  fprintf('%15.2f', Y(i, 1));

  time = tic;
  kt.solve(PP, zeros(kt.sdim, 1));
  Y(i, 2) = toc(time) * monteCarloSamples;

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

legend('Proposed Framework', 'One Monte Carlo');
