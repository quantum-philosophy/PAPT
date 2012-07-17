init;

hsConfig = Utils.resolvePath('hotspot.config');
floorplan = Utils.resolvePath('dummy2.flp');
powerTrace = Utils.resolvePath('dummy2.ptrace');

hsLine = 'sampling_intvl 1e-3';

ks = KuttaSpot(floorplan, hsConfig, hsLine);

fprintf('Sampling interval:   %.2e s\n', ks.dt);
fprintf('Ambient temperature: %.2f K\n', ks.Tamb);
fprintf('Number of nodes:     %d\n', ks.nodes);
fprintf('Number of cores:     %d\n', ks.cores);

Pdyn = dlmread(powerTrace, '', 1, 0)';

steps = size(Pdyn, 2);

fprintf('Number of steps:     %d\n', steps);
fprintf('Total time:          %.2f s\n', ks.dt * steps);

t = tic;
T = ks.solve(Pdyn, zeros(ks.sdim, 1)) + Constants.zeroKelvin;
fprintf('Simulation time:     %.2f s\n', toc(t));

time = (1:steps) * ks.dt;

figure;
for i = 1:ks.cores
  color = Utils.pickColor(i);
  line(time, T(i, :), 'Color', color);
end
