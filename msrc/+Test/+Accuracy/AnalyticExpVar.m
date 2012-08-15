init;

kc = Config('steps', 100, 'assessmentMethod', 'Kutta');
kt = Test.constructMonteCarlo(kc);
[ kExp, kVar, kRaw ] = Test.sampleMonteCarlo(kt, kc);

ac = Config('steps', 100, 'assessmentMethod', 'Analytic');
an = Test.constructMonteCarlo(ac);
[ aExp, aVar, aRaw ] = Test.sampleMonteCarlo(an, ac);

time = kc.timeLine;

%% Temperature analysis with Kutta
%

kf = Stats.drawEvolution(time, kExp, kVar);
title(sprintf('%d-sample Monte Carlo with Kutta', kc.monteCarloSamples));

%% Temperature analysis with Analytic
%

af = Stats.drawEvolution(time, aExp, aVar);
title(sprintf('%d-sample Monte Carlo with Analytic', ac.monteCarloSamples));

Utils.evenScale(kf, af);

%% Comparison
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

ExpNRMSE = Stats.NRMSE(kExp, aExp) * 100;
VarNRMSE = Stats.NRMSE(kVar, aVar) * 100;

fprintf('NRMSE(Exp): %.2f%%\n', ExpNRMSE);
fprintf('NRMSE(Var): %.2f%%\n', VarNRMSE);
