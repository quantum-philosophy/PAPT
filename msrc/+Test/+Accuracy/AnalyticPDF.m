init;

kc = Config('assessmentMethod', 'Kutta');
kt = Test.constructMonteCarlo(kc);
[ ~, ~, kRaw ] = Test.sampleMonteCarlo(kt, kc);

ac = Config('assessmentMethod', 'Analytic');
an = Test.constructMonteCarlo(ac);
[ ~, ~, aRaw ] = Test.sampleMonteCarlo(an, ac);

[ globalError, localError ] = ...
  Utils.comparePDF(kRaw, aRaw, kc.timeLine, { 'Kutta', 'Analytic' });

fprintf('\n');
fprintf('%10s%20s\n', 'Core', 'Average NRMSE, %');

for i = 1:kc.cores
  fprintf('%10d%20.2f\n', i, mean(localError(i, :)) * 100);
end

fprintf('\n');
fprintf('Global error: %.2f%%\n', globalError * 100);
