clear all;
clc;

floorplan = Utils.resolvePath('dummy.flp');
hsConfig = Utils.resolvePath('hotspot.config');
hsLine = 'sampling_intvl 1e-3';

hotspot = HotSpot(floorplan, hsConfig, hsLine);

fprintf('Sampling interval:   %.2e s\n', hotspot.dt);
fprintf('Ambient temperature: %.2f K\n', hotspot.Tamb);
fprintf('Number of nodes:     %d\n', hotspot.nodes);
fprintf('Number of cores:     %d\n', hotspot.cores);

powerTrace = Utils.resolvePath('dummy.ptrace');
P = dlmread(powerTrace, '', 1, 0)';

steps = size(P, 2);

fprintf('Number of steps:     %d\n', steps);
fprintf('Total time:          %.2f s\n', hotspot.dt * steps);

T = hotspot.solve(P) + Constants.zeroKelvin;

time = (1:steps) * hotspot.dt;

plot(time, T);
