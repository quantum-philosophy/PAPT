init;

[ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(2);

hs = HotSpot.Analytic(floorplan, config, configLine);

fprintf('Sampling interval:   %.2e s\n', hs.dt);
fprintf('Ambient temperature: %.2f K\n', hs.Tamb);
fprintf('Number of nodes:     %d\n', hs.nodes);
fprintf('Number of cores:     %d\n', hs.cores);

Pdyn = dlmread(powerTrace, '', 1, 0)';

steps = size(Pdyn, 2);

fprintf('Number of steps:     %d\n', steps);
fprintf('Total time:          %.2f s\n', hs.dt * steps);

t = tic;
[ ExpT, VarT ] = MonteCarlo.perform3D( ...
  @(rvs) hs.solve(Pdyn, rvs), [ hs.sdim, hs.cores, steps ]);
t = toc(t);
fprintf('Simulation time:     %.2f s\n', t);

ExpT = ExpT + Constants.zeroKelvin;
StdT = sqrt(VarT);

time = (1:steps) * hs.dt;

figure;
for i = 1:hs.cores
  color = Utils.pickColor(i);
  line(time, ExpT(i, :), 'Color', color);
  line(time, ExpT(i, :) + StdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, ExpT(i, :) - StdT(i, :), 'Color', color, 'LineStyle', '--');
end
title(sprintf('Monte Carlo (%.2f s)', t));
