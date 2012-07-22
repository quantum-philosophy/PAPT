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
T = hs.solve(Pdyn, zeros(hs.sdim, 1)) + Constants.zeroKelvin;
t = toc(t);
fprintf('Simulation time:     %.2f s\n', t);

time = (1:steps) * hs.dt;

figure;
title(sprintf('HotSpot with Analytical Solution (%.2f s)', t));
for i = 1:hs.cores
  color = Utils.pickColor(i);
  line(time, T(i, :), 'Color', color);
end
