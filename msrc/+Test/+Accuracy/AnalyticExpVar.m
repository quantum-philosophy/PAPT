init;

kc = Config('assessmentMethod', 'Kutta', 'steps', 100);
kt = Test.constructMonteCarlo(kc);
[ kExp, kVar, kRaw ] = Test.sampleMonteCarlo(kt, kc);

ac = Config('assessmentMethod', 'Analytic', 'steps', 100);
an = Test.constructMonteCarlo(ac);
[ aExp, aVar, aRaw ] = Test.sampleMonteCarlo(an, ac);

time = kc.timeLine;

%% Temperature analysis with Kutta.
%

Utils.plotExpStd(time, kExp, kVar);
title(sprintf('%d-sample Monte Carlo with Kutta', kc.monteCarloSamples));

%% Temperature analysis with Analytic.
%

Utils.plotExpStd(time, aExp, aVar);
title(sprintf('%d-sample Monte Carlo with Analytic', ac.monteCarloSamples));

%% Comparison of the methods.
%

figure;
subplot(2, 1, 1);
title('Difference of Expectations');
for i = 1:kc.cores
  color = Utils.pickColor(i);
  line(time, kExp(i, :) - aExp(i, :), 'Color', color);
end
ylabel('Exp(Kutta) - Exp(Analytic), C');
xlabel('Time, s');

subplot(2, 1, 2);
title('Difference of Variances');
for i = 1:kc.cores
  color = Utils.pickColor(i);
  line(time, kVar(i, :) - aVar(i, :), 'Color', color);
end
ylabel('Var(Kutta) - Var(Analytic), C^2');
xlabel('Time, s');

ExpNRMSE = Utils.NRMSE(kExp, aExp) * 100;
VarNRMSE = Utils.NRMSE(kVar, aVar) * 100;

fprintf('NRMSE(Exp): %.2f%%\n', ExpNRMSE);
fprintf('NRMSE(Var): %.2f%%\n', VarNRMSE);
