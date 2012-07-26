init;

cores = 2;
steps = 100;
divide = 2;
samples = 2e4;
method = 'Kutta';
dt = 1e-3;

X = [ 1 2 3 4 5 6 7 8 9 10 ];

pick = length(X);

[ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(cores);
configLine = sprintf('sampling_intvl %.2e', dt);

fprintf('Number of cores: %d\n', cores);
fprintf('Number of steps: %d\n', steps);
fprintf('Sampling interval: %.2e s\n', dt);

Pdyn = Utils.replicate(dlmread(powerTrace, '', 1, 0)', steps);

%% Monte Carlo sampling.
%

hs = HotSpot.(method)(floorplan, config, configLine);

assert(hs.dt == dt, 'The sampling interval is invalid.');

filename = sprintf('MonteCarlo_%s_nc%d_ns%d_f%d_mcs%d.mat', ...
  method, cores, steps, 1 / dt, samples);
filename = Utils.resolvePath(filename, 'cache');

if exist(filename)
  load(filename);
else
  t = tic;
  [ mcExpT, mcVarT ] = MonteCarlo.perform3D( ...
    @(rvs) hs.solve(Pdyn, rvs), [ hs.sdim, cores, steps ], samples);
  t = toc(t);

  save(filename, 't', 'mcExpT', 'mcVarT');
end

fprintf('MC simulation time: %.2f s\n', t);

mcExpT = Utils.toCelsius(mcExpT);
mcStdT = sqrt(mcVarT);

tmc = t;

%% Polynomial Chaos expansion.
%

Pdyn0 = Pdyn;
Pdyn = zeros(cores, steps * divide);

I = ((1:steps) - 1) * divide + 1;
for i = 1:divide
  Pdyn(:, I + i - 1) = Pdyn0;
end
I = I + (divide - 1);

configLine = sprintf('sampling_intvl %.2e', dt / divide);

fprintf('Extended number of steps: %d\n', steps * divide);
fprintf('Modified sampling interval: %.2e s\n', dt / divide);

fprintf('%15s%15s%15s%15s\n', 'Order', 'Time, s', 'NRMSE(Exp), %', 'NRMSE(Var), %');

count = length(X);

ExpNRMSE = zeros(count, 1);
VarNRMSE = zeros(count, 1);

for i = 1:count
  order = X(i);

  hs = HotSpot.Chaos(floorplan, config, configLine, order);
  assert(hs.dt == dt / divide, 'The sampling interval is invalid.');

  t = tic;
  [ ExpT, VarT ] = hs.solve(Pdyn);
  t = toc(t);

  ExpT = ExpT(:, I);
  VarT = VarT(:, I);

  ExpT = Utils.toCelsius(ExpT);
  StdT = sqrt(VarT);

  ExpNRMSE(i) = Utils.NRMSE(mcExpT, ExpT, 5) * 100;
  VarNRMSE(i) = Utils.NRMSE(mcVarT, VarT, 5) * 100;

  fprintf('%15d%15.2f%15.2f%15.2f\n', order, t, ...
    ExpNRMSE(i), VarNRMSE(i));

  if pick == i
    pcExpT = ExpT;
    pcStdT = StdT;
    pcVarT = VarT;
  end
end

time = (1:steps) * dt;

mcf = figure;
title(sprintf('%d-sample Monte Carlo (%.2f s)', samples, tmc));
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, mcExpT(i, :), 'Color', color);
  line(time, mcExpT(i, :) + mcStdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, mcExpT(i, :) - mcStdT(i, :), 'Color', color, 'LineStyle', '--');
end

pcf = figure;
title(sprintf('%d-order Polynomial Chaos (%.2f s)', X(pick), t));
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, pcExpT(i, :), 'Color', color);
  line(time, pcExpT(i, :) + pcStdT(i, :), 'Color', color, 'LineStyle', '--');
  line(time, pcExpT(i, :) - pcStdT(i, :), 'Color', color, 'LineStyle', '--');
end

%% Make the plots match each other.
%
lim1 = ylim;
figure(mcf);
lim2 = ylim;
lim = [ min(lim1(1), lim2(1)) max(lim1(2), lim2(2)) ];
ylim(lim);
figure(pcf);
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

figure;
subplot(2, 1, 1);
title('Convergence of Expectations');
color = Utils.pickColor(1);
line(X, ExpNRMSE, 'Color', color);
xlabel('Order');
ylabel('NRMSE(Exp), %');

subplot(2, 1, 2);
title('Convergence of Variances');
color = Utils.pickColor(1);
line(X, VarNRMSE, 'Color', color);
xlabel('Order');
ylabel('NRMSE(Var), %');
