init;

floorplan = Utils.resolvePath('dummy.flp');
hsConfig = Utils.resolvePath('hotspot.config');
hsLine = 'sampling_intvl 1e-3';

hs = HotSpot(floorplan, hsConfig, hsLine);

fprintf('Sampling interval:   %.2e s\n', hs.dt);
fprintf('Ambient temperature: %.2f K\n', hs.Tamb);
fprintf('Number of nodes:     %d\n', hs.nodes);
fprintf('Number of cores:     %d\n', hs.cores);

powerTrace = Utils.resolvePath('dummy.ptrace');
Pdyn = dlmread(powerTrace, '', 1, 0)';

steps = size(Pdyn, 2);

fprintf('Number of steps:     %d\n', steps);
fprintf('Total time:          %.2f s\n', hs.dt * steps);

leakage = Leakage.constructBasedOnDynamic(Pdyn);

mc = MonteCarlo(hs, Pdyn, leakage);

t = tic;
[ ExpT, VarT ] = MonteCarlo.perform3D(...
  @mc.evaluate, [ hs.cores, hs.cores, steps ]);
t = toc(t);
fprintf('Simulation time:     %.2f s\n', t);

ExpT = ExpT + Constants.zeroKelvin;
StdT = VarT.^(1/2);

time = (1:steps) * hs.dt;

figure;
for i = 1:hs.cores
  color = Utils.pickColor(i);
  line(time, ExpT(i, :), 'Color', color);
  line(time, ExpT(i, :) + StdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, ExpT(i, :) - StdT(i, :), 'Color', color, 'LineStyle', '--');
end
title(sprintf('Monte-Carlo (%.2f s)', t));
