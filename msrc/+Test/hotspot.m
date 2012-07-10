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

t = tic;
T = hs.solve(Pdyn, leakage) + Constants.zeroKelvin;
fprintf('Simulation time:     %.2f s\n', toc(t));

time = (1:steps) * hs.dt;

plot(time, T);
