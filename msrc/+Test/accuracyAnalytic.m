init;

cores = 4;
steps = 100;
samples = 1e4;

[ floorplan, powerTrace, config, configLine ] = Utils.resolveTest(cores);

fprintf('Number of cores: %d\n', cores);
fprintf('Number of steps: %d\n', steps);

Pdyn = Utils.replicate(dlmread(powerTrace, '', 1, 0)', steps);

%% Monte Carlo sampling.
%

methodSet = { 'Analytic', 'Kutta' };

ExpT = {};
VarT = {};
Time = {};

for i = 1:length(methodSet)
  method = methodSet{i};

  hs = HotSpot.(method)(floorplan, config, configLine);

  stamp = sprintf('%s_nc%d_ns%d_f%d', method, cores, steps, 1 / hs.dt);

  [ mcExpT, mcVarT, t ] = MonteCarlo.perform3D( ...
    @(rvs) hs.solve(Pdyn, rvs), [ hs.sdim, cores, steps ], samples, stamp);

  fprintf('MC with %s simulation time:  %.2f s\n', method, t);

  mcExpT = Utils.toCelsius(mcExpT);
  mcStdT = sqrt(mcVarT);

  ExpT{end + 1} = mcExpT;
  VarT{end + 1} = mcVarT;
  Time{end + 1} = t;

  time = (1:steps) * hs.dt;

  figure;
  title(sprintf('%d-sample Monte Carlo with %s (%.2f s)', samples, method, t));
  for i = 1:cores
    color = Utils.pickColor(i);
    line(time, mcExpT(i, :), 'Color', color);
    line(time, mcExpT(i, :) + mcStdT(i, :), 'Color', color, 'LineStyle', '--');
    line(time, mcExpT(i, :) - mcStdT(i, :), 'Color', color, 'LineStyle', '--');
  end
end

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Difference of Expectations');
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, ExpT{2}(i, :) - ExpT{1}(i, :), 'Color', color);
end

subplot(2, 1, 2);
title('Difference of Variances');
for i = 1:cores
  color = Utils.pickColor(i);
  line(time, VarT{2}(i, :) - VarT{1}(i, :), 'Color', color);
end

ExpNRMSE = Utils.NRMSE(ExpT{2}, ExpT{1}, 5) * 100;
VarNRMSE = Utils.NRMSE(VarT{2}, VarT{1}, 5) * 100;

fprintf('NRMSE(Exp): %.2f %%\n', ExpNRMSE);
fprintf('NRMSE(Var): %.2f %%\n', VarNRMSE);
