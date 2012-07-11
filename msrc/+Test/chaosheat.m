init;

hsConfig = Utils.resolvePath('hotspot.config');
floorplan = Utils.resolvePath('dummy2.flp');
powerTrace = Utils.resolvePath('dummy2.ptrace');

hsLine = 'sampling_intvl 1e-3';

ch = ChaosHeat(floorplan, hsConfig, hsLine);

fprintf('Sampling interval:   %.2e s\n', ch.dt);
fprintf('Ambient temperature: %.2f K\n', ch.Tamb);
fprintf('Number of nodes:     %d\n', ch.nodes);
fprintf('Number of cores:     %d\n', ch.cores);

Pdyn = dlmread(powerTrace, '', 1, 0)';

steps = size(Pdyn, 2);

fprintf('Number of steps:     %d\n', steps);
fprintf('Total time:          %.2f s\n', ch.dt * steps);

t = tic;
[ ExpT, VarT ] = ch.solve(Pdyn);
t = toc(t);
fprintf('Simulation time:     %.2f s\n', t);

ExpT = ExpT + Constants.zeroKelvin;
StdT = VarT.^(1/2);

time = (1:steps) * ch.dt;

figure;
for i = 1:ch.cores
  color = Utils.pickColor(i);
  line(time, ExpT(i, :), 'Color', color);
  line(time, ExpT(i, :) + StdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, ExpT(i, :) - StdT(i, :), 'Color', color, 'LineStyle', '--');
end
title(sprintf('Polynomial Chaos (%.2f s)', t));
