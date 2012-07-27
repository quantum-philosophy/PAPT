init;

c = Test.config();
c.adjustPowerSteps(100);
c.display();

time = c.timeLine;

%% Temperature analysis with Kutta.
%

[ kExp, kVar, kRaw, kTime ] = Test.computeKutta(c);
fprintf('MC with Kutta simulation time: %.2f s\n', kTime);

std = sqrt(kVar);

figure;
title(sprintf('%d-sample Monte Carlo with Kutta (%.2f s)', c.samples, kTime));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, kExp(i, :), 'Color', color);
  line(time, kExp(i, :) + std(i, :), 'Color', color, 'LineStyle', '--');
  line(time, kExp(i, :) - std(i, :), 'Color', color, 'LineStyle', '--');
end

%% Temperature analysis with Analytic.
%

[ aExp, aVar, aRaw, aTime ] = Test.computeAnalytic(c);
fprintf('MC with Analytic simulation time: %.2f s\n', aTime);

std = sqrt(aVar);

figure;
title(sprintf('%d-sample Monte Carlo with Analytic (%.2f s)', c.samples, aTime));
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, aExp(i, :), 'Color', color);
  line(time, aExp(i, :) + std(i, :), 'Color', color, 'LineStyle', '--');
  line(time, aExp(i, :) - std(i, :), 'Color', color, 'LineStyle', '--');
end

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Difference of Expectations');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, kExp(i, :) - aExp(i, :), 'Color', color);
end

subplot(2, 1, 2);
title('Difference of Variances');
for i = 1:c.cores
  color = Utils.pickColor(i);
  line(time, kVar(i, :) - aVar(i, :), 'Color', color);
end

ExpNRMSE = Utils.NRMSE(kExp, aExp) * 100;
VarNRMSE = Utils.NRMSE(kVar, aVar) * 100;

fprintf('NRMSE(Exp): %.2f%%\n', ExpNRMSE);
fprintf('NRMSE(Var): %.2f%%\n', VarNRMSE);
