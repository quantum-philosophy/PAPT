init;

c = Test.config();
c.adjustPowerSteps(100);
c.display();

time = c.timeLine;

[ kExp, kVar, kRaw, kTime ] = Test.computeKutta(c);
[ aExp, aVar, aRaw, aTime ] = Test.computeAnalytic(c);

error = Utils.comparePDF(c.timeLine, kRaw, aRaw, { 'Kutta', 'Analytic' });

fprintf('%10s%20s\n', 'Core', 'Average NRMSE');
for i = 1:c.cores
  fprintf('%10d%20.2e\n', i, mean(error(i, :)));
end
