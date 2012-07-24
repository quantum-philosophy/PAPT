init;

cores = 2;
steps = 100;
order = 3;
samples = 1e4;

[ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(cores);

fprintf('Number of cores:     %d\n', cores);
fprintf('Number of steps:     %d\n', steps);
fprintf('PC order:            %d\n', order);

Pdyn = Utils.replicate(dlmread(powerTrace, '', 1, 0)', steps);

%% Polynomial Chaos expansion.
%

hs = HotSpot.Chaos(floorplan, config, configLine, order);

t = tic;
[ pcExpT, pcVarT ] = hs.solve(Pdyn);
t = toc(t);
fprintf('PC simulation time:  %.2f s\n', t);

pcExpT = Utils.toCelsius(pcExpT);
pcStdT = sqrt(pcVarT);

time = (1:steps) * hs.dt;

pcf = figure;
title(sprintf('%d-order Polynomial Chaos (%.2f s)', order, t));
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, pcExpT(i, :), 'Color', color);
  line(time, pcExpT(i, :) + pcStdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, pcExpT(i, :) - pcStdT(i, :), 'Color', color, 'LineStyle', '--');
end

%% Monte Carlo sampling.
%

filename = sprintf('MonteCarlo_c%d_s%d_s%d.mat', cores, steps, samples);
filename = Utils.resolvePath(filename, 'cache');

if exist(filename)
  load(filename);
else
  hs = HotSpot.Kutta(floorplan, config, configLine);

  t = tic;
  [ mcExpT, mcVarT ] = MonteCarlo.perform3D( ...
    @(rvs) hs.solve(Pdyn, rvs), [ hs.sdim, cores, steps ], samples);
  t = toc(t);

  save(filename, 't', 'mcExpT', 'mcVarT');
end

fprintf('MC simulation time:  %.2f s\n', t);

mcExpT = Utils.toCelsius(mcExpT);
mcStdT = sqrt(mcVarT);

mcf = figure;
title(sprintf('%d-sample Monte Carlo (%.2f s)', samples, t));
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, mcExpT(i, :), 'Color', color);
  line(time, mcExpT(i, :) + mcStdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, mcExpT(i, :) - mcStdT(i, :), 'Color', color, 'LineStyle', '--');
end

%% Make the plots match each other.
%
mclim = ylim;
figure(pcf);
pclim = ylim;
lim = [ min(pclim(1), mclim(1)) max(pclim(2), mclim(2)) ];
ylim(lim);
figure(mcf);
ylim(lim);

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Difference of Expectations');
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, mcExpT(i, :) - pcExpT(i, :), 'Color', color);
end

subplot(2, 1, 2);
title('Difference of Variances');
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, mcVarT(i, :) - pcVarT(i, :), 'Color', color);
end
