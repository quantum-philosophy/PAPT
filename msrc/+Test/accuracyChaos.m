init;

cores = 4;
steps = 100;
samples = 1e4;
method = 'Kutta';

X = [ 1 2 3 4 5 6 ];

pick = length(X);

[ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(cores);

fprintf('Number of cores: %d\n', cores);
fprintf('Number of steps: %d\n', steps);
fprintf('Number of samples: %d\n', samples);

Pdyn = Utils.replicate(dlmread(powerTrace, '', 1, 0)', steps);

%% Monte Carlo sampling.
%

hs = HotSpot.(method)(floorplan, config, configLine);

stamp = sprintf('%s_nc%d_ns%d_f%d', method, cores, steps, 1 / hs.dt);

[ mcExpT, mcVarT, t ] = MonteCarlo.perform3D( ...
  @(rvs) hs.solve(Pdyn, rvs), [ hs.sdim, cores, steps ], samples, stamp);

fprintf('MC simulation time: %.2f s\n', t);

mcExpT = Utils.toCelsius(mcExpT);
mcStdT = sqrt(mcVarT);

tmc = t;

%% Polynomial Chaos expansion.
%

fprintf('%15s%15s%15s%15s\n', 'Order', 'Time, s', 'NRMSE(Exp), %', 'NRMSE(Var), %');

count = length(X);

ExpNRMSE = zeros(count, 1);
VarNRMSE = zeros(count, 1);

for i = 1:count
  order = X(i);

  hs = HotSpot.Chaos(floorplan, config, configLine, order);

  t = tic;
  [ ExpT, VarT ] = hs.solve(Pdyn);
  t = toc(t);

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

time = (1:steps) * hs.dt;

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
