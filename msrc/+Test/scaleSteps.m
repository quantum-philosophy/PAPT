init;

floorplan = Utils.resolvePath('dummy2.flp');
config = Utils.resolvePath('hotspot.config');
powerTrace = Utils.resolvePath('dummy2.ptrace');

configLine = 'sampling_intvl 1e-3';

monteCarloSamples = 1e4;

%% Initialize the solver.
%
ch = HotSpot.Chaos(floorplan, config, configLine);
kt = HotSpot.Kutta(floorplan, config, configLine);

%% Define the needed measurements.
%
X = [ 10, 100, 1000, 10000, 100000 ];

%% Read the specified template.
%
P = dlmread(powerTrace, '', 1, 0)';
[ cores, temps ] = size(P);

fprintf('%15s%15s%15s%15s\n', 'Steps', 'Chaos, s', 'Kutta, s', 'Speedup, x');

Y = zeros(length(X), 2);

for i = 1:length(X)
  steps = X(i);

  fprintf('%15d', steps);

  %% Construct a power profile of specific size.
  %
  PP = zeros(cores, steps);
  packed = 0;
  while (packed < steps)
    topack = min(steps - packed, temps);
    PP(:, (packed + 1):(packed + topack)) = P(:, 1:topack);
    packed = packed + topack;
  end

  %% Perform the analysis.
  %
  time = tic;
  ch.solve(PP);
  Y(i, 1) = toc(time);

  fprintf('%15.2f', Y(i, 1));

  time = tic;
  kt.solve(PP, zeros(kt.sdim, 1));
  Y(i, 2) = toc(time) * monteCarloSamples;

  fprintf('%15.2f', Y(i, 2));

  fprintf('%15.2e', Y(i, 2) / Y(i, 1));
  fprintf('\n');
end

figure;
line(X, Y(:, 1), 'Color', Utils.pickColor(1), 'Marker', 'o')
line(X, Y(:, 2), 'Color', Utils.pickColor(2), 'Marker', 'x')
title('Scaling with Number of Steps');
xlabel('Steps');
ylabel('Time, s');

legend('Proposed Framework', 'One Monte Carlo');
